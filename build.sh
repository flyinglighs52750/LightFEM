#!/usr/bin/env bash
# 构建 LightFEM（Release）。需要 CMake >= 3.16 与支持 C++17 的编译器。
#
# 用法：
#   ./build.sh                         # 配置并构建到 ./build
#   BUILD_DIR=build-release ./build.sh # 指定构建目录
#   CMAKE_GENERATOR=Ninja ./build.sh   # 指定 CMake 生成器
#
# Windows 下会先用 vswhere 定位 Visual Studio 并直接注入 MSVC 的 INCLUDE/LIB
# （不依赖 vcvars64.bat，也不需要以管理员身份运行），再用 Ninja 构建。
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT/build}"

find_cmake() {
    if command -v cmake >/dev/null 2>&1; then
        command -v cmake
        return 0
    fi
    local c
    for c in \
        "/c/Program Files/Microsoft Visual Studio/2022/Enterprise/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe" \
        "/c/Program Files/Microsoft Visual Studio/2022/Professional/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe" \
        "/c/Program Files/Microsoft Visual Studio/2022/Community/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe" \
        "/d/VS/VS2022/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe"; do
        if [ -x "$c" ]; then
            printf '%s\n' "$c"
            return 0
        fi
    done
    printf '未找到 cmake，请安装 CMake 3.16 或更高版本，或将其加入 PATH\n' >&2
    return 1
}

# 用 vswhere 定位 Visual Studio，并导出 MSVC 编译所需的环境变量。
prepare_msvc() {
    local vswhere="/c/Program Files (x86)/Microsoft Visual Studio/Installer/vswhere.exe"
    [ -x "$vswhere" ] || return 1

    local vs_path
    vs_path="$("$vswhere" -latest -products '*' \
        -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 \
        -property installationPath 2>/dev/null | tr -d '\r' | head -1)"
    [ -n "$vs_path" ] || return 1

    local win_vs
    win_vs="$(printf '%s' "$vs_path" | sed 's|^\([A-Za-z]\):|/\L\1|; s|\\|/|g')"

    local msvc_dir
    msvc_dir="$(ls -d "$win_vs"/VC/Tools/MSVC/*/ 2>/dev/null | sort -V | tail -1)"
    [ -n "$msvc_dir" ] || return 1
    msvc_dir="${msvc_dir%/}"

    local sdk_root="/c/Program Files (x86)/Windows Kits/10"
    local sdk_ver
    sdk_ver="$(ls -d "$sdk_root"/Include/*/ 2>/dev/null | sort -V | tail -1)"
    [ -n "$sdk_ver" ] || return 1
    sdk_ver="$(basename "$sdk_ver")"

    # MSVC 的 cl.exe / link.exe 只认 Windows 风格路径，这里统一转换。
    local msvc_win sdk_win
    msvc_win="$(cygpath -m "$msvc_dir")"
    sdk_win="$(cygpath -m "$sdk_root")"

    export INCLUDE="$msvc_win/include;$sdk_win/Include/$sdk_ver/ucrt;$sdk_win/Include/$sdk_ver/um;$sdk_win/Include/$sdk_ver/shared;$sdk_win/Include/$sdk_ver/winrt"
    export LIB="$msvc_win/lib/x64;$sdk_win/Lib/$sdk_ver/ucrt/x64;$sdk_win/Lib/$sdk_ver/um/x64"
    export PATH="$msvc_dir/bin/Hostx64/x64:$sdk_root/bin/$sdk_ver/x64:$win_vs/Common7/IDE/CommonExtensions/Microsoft/CMake/Ninja:$PATH"
    MSVC_CL="$msvc_win/bin/Hostx64/x64/cl.exe"
    export MSVC_CL
    return 0
}

CMAKE="$(find_cmake)"

case "$(uname -s)" in
    MINGW* | MSYS* | CYGWIN*)
        if [ -n "${CMAKE_GENERATOR:-}" ] && [ "${CMAKE_GENERATOR}" != "Ninja" ]; then
            "$CMAKE" -S "$ROOT" -B "$BUILD_DIR" -G "$CMAKE_GENERATOR" -A x64
            "$CMAKE" --build "$BUILD_DIR" --config Release --parallel
        elif prepare_msvc; then
            "$CMAKE" -S "$ROOT" -B "$BUILD_DIR" -G Ninja \
                -DCMAKE_BUILD_TYPE=Release \
                -DCMAKE_C_COMPILER="$(dirname "$MSVC_CL")/cl.exe" \
                -DCMAKE_CXX_COMPILER="$MSVC_CL"
            "$CMAKE" --build "$BUILD_DIR"
        else
            printf '未找到可用的 Visual Studio C++ 工具集，改用默认生成器\n' >&2
            "$CMAKE" -S "$ROOT" -B "$BUILD_DIR"
            "$CMAKE" --build "$BUILD_DIR" --config Release --parallel
        fi
        ;;
    *)
        "$CMAKE" -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
        "$CMAKE" --build "$BUILD_DIR" --parallel
        ;;
esac

printf '\n构建完成：%s\n' "$BUILD_DIR"
