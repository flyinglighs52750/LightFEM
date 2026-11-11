#!/bin/bash
# 构建 Release 版内核（用于发行包里的 lightfem.exe）。
# 与 build.sh 同样的手工拼 INCLUDE / LIB，只是目标目录换成 build-release。
set -e
MSVC_DIR="D:/VS/VS2022/VC/Tools/MSVC/14.44.35207"
SDK_VER="10.0.26100.0"
SDK="C:/Program Files (x86)/Windows Kits/10"
export INCLUDE="$MSVC_DIR/include;$SDK/Include/$SDK_VER/ucrt;$SDK/Include/$SDK_VER/um;$SDK/Include/$SDK_VER/shared;$SDK/Include/$SDK_VER/winrt"
export LIB="$MSVC_DIR/lib/x64;$SDK/Lib/$SDK_VER/ucrt/x64;$SDK/Lib/$SDK_VER/um/x64"
NINJA="D:/VS/VS2022/Common7/IDE/CommonExtensions/Microsoft/CMake/Ninja/ninja.exe"
cd "$(dirname "$0")/build-release"
"$NINJA" "$@"
