@echo off
rem 构建 LightFEM（Release）。需要 CMake >= 3.16 与支持 C++17 的编译器。
rem 用法：build.bat  [可选：附加的 cmake --build 参数]
rem 建议在 “x64 Native Tools Command Prompt for VS 2022” 中运行；若 cmake 已在
rem PATH 中，则任意命令行均可。

setlocal
set "ROOT=%~dp0"
if "%BUILD_DIR%"=="" set "BUILD_DIR=%ROOT%build"

where cmake >nul 2>&1
if not errorlevel 1 (
  set "CMAKE=cmake"
) else if exist "%ProgramFiles%\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" (
  set "CMAKE=%ProgramFiles%\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
) else if exist "%ProgramFiles%\Microsoft Visual Studio\2022\Professional\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" (
  set "CMAKE=%ProgramFiles%\Microsoft Visual Studio\2022\Professional\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
) else if exist "D:\VS\VS2022\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" (
  set "CMAKE=D:\VS\VS2022\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
) else (
  echo 未找到 cmake，请安装 CMake 3.16 或更高版本，或将其加入 PATH
  exit /b 1
)

"%CMAKE%" -S "%ROOT%" -B "%BUILD_DIR%" -G "Visual Studio 17 2022" -A x64
if errorlevel 1 exit /b 1

"%CMAKE%" --build "%BUILD_DIR%" --config Release --parallel %*
if errorlevel 1 exit /b 1

echo.
echo 构建完成，可执行文件位于 %BUILD_DIR%\Release
endlocal
