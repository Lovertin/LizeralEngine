# LizeralEngine C++ Project

A modern C++ project with CMake support.

## Project Structure

```
.
├── CMakeLists.txt      # CMake configuration
├── src/               # Source files
│   ├── main.cpp
│   └── example.cpp
├── include/           # Header files
│   └── example.h
├── build/             # Build directory (generated)
├── build.bat          # Windows build script
└── build.sh           # Unix/Linux/macOS build script
```

## Prerequisites

- CMake 3.10 or higher
- C++17 compatible compiler
  - Windows: Visual Studio 2022 or MinGW
  - Linux: GCC 7+ or Clang 5+
  - macOS: Xcode 10+

## Building the Project

### Windows
```bash
build.bat
```

Or manually:
```bash
mkdir build
cd build
cmake .. -G "Visual Studio 17 2022" -A x64
cmake --build . --config Release
```

### Unix/Linux/macOS
```bash
chmod +x build.sh
./build.sh
```

Or manually:
```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

## Running the Executable

After building, the executable will be located at:
- Windows: `build/bin/Release/LizeralEngine.exe`
- Unix-like: `build/bin/LizeralEngine`

Run it with:
```bash
# Windows
build\bin\Release\LizeralEngine.exe

# Unix-like
./build/bin/LizeralEngine
```

## Adding New Files

1. Add source files to `src/` directory
2. Add header files to `include/` directory
3. CMake will automatically detect `.cpp`, `.cxx`, `.cc` files in `src/`
4. CMake will automatically detect `.h`, `.hpp`, `.hxx` files in `include/`

## CMake Configuration

Key features of the CMake configuration:
- C++17 standard required
- Strict compiler warnings (treated as errors)
- Organized output directories (bin/, lib/)
- Cross-platform compiler options
- Install targets

## Customization

Edit `CMakeLists.txt` to:
- Change project name or version
- Add dependencies with `find_package()`
- Add subdirectories with `add_subdirectory()`
- Modify compiler flags
- Add libraries with `target_link_libraries()`

## License

This is a template project. Feel free to modify as needed.
