#!/usr/bin/env bash
# 构建 Release 版内核（用于发行包）。等价于 BUILD_DIR=build-release ./build.sh。
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export BUILD_DIR="$ROOT/build-release"
exec "$ROOT/build.sh" "$@"
