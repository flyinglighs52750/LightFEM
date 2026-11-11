@echo off
call "D:\VS\VS2022\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 (
  echo vcvars64.bat failed
  exit /b 1
)
cd /d "%~dp0build-ninja"
"D:\VS\VS2022\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja\ninja.exe" %*
