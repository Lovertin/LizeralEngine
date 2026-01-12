#!/bin/bash
# Build script for Unix-like systems (Linux, macOS)

echo "Setting up build directory..."
mkdir -p build
cd build

echo "Configuring project with CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release

if [ $? -ne 0 ]; then
    echo "CMake configuration failed!"
    exit 1
fi

echo "Building project..."
cmake --build . --config Release

if [ $? -ne 0 ]; then
    echo "Build failed!"
    exit 1
fi

echo "Build successful!"
echo ""
echo "Executable location: build/bin/LizeralEngine"
echo ""
cd ..
