@echo off
REM Build script for Windows

echo Setting up build directory...
if not exist "build" mkdir build
cd build

echo Configuring project with CMake...
cmake .. -G "Visual Studio 17 2022" -A x64

if %ERRORLEVEL% NEQ 0 (
    echo CMake configuration failed!
    exit /b 1
)

echo Building project...
cmake --build . --config Release

if %ERRORLEVEL% NEQ 0 (
    echo Build failed!
    exit /b 1
)

echo Build successful!
echo.
echo Executable location: build\bin\Release\LizeralEngine.exe
echo.
cd ..
