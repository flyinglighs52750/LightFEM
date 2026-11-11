#!/bin/bash
# 用 MSVC 环境构建 LightFEM（ninja）。
# 不调用 vcvars64.bat：它内部会调 reg.exe 探测 SDK，而 reg 被安全策略拦。
# 这里直接按已知路径拼 INCLUDE / LIB。
set -e
MSVC_DIR="D:/VS/VS2022/VC/Tools/MSVC/14.44.35207"
SDK_VER="10.0.26100.0"
SDK="C:/Program Files (x86)/Windows Kits/10"
export INCLUDE="$MSVC_DIR/include;$SDK/Include/$SDK_VER/ucrt;$SDK/Include/$SDK_VER/um;$SDK/Include/$SDK_VER/shared;$SDK/Include/$SDK_VER/winrt"
export LIB="$MSVC_DIR/lib/x64;$SDK/Lib/$SDK_VER/ucrt/x64;$SDK/Lib/$SDK_VER/um/x64"
NINJA="D:/VS/VS2022/Common7/IDE/CommonExtensions/Microsoft/CMake/Ninja/ninja.exe"
cd "$(dirname "$0")/build-ninja"
"$NINJA" "$@"
