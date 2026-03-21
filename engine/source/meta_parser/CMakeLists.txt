cmake_minimum_required(VERSION 3.10)
project(LizeralParser)

set(TARGET_NAME LizeralParser)

# 开启文件夹支持（用于IDE分类）
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# =========================================================
# 1. 源文件收集 (Source Collection)
# =========================================================
file(GLOB_RECURSE HEADERS "*.h")
file(GLOB_RECURSE ALL_SOURCES "*.cpp")

# 过滤掉CMake生成的临时检测文件
set(SOURCES "")
foreach(source ${ALL_SOURCES})
    if(NOT source MATCHES "CMakeFiles.*CompilerIdCXX")
        list(APPEND SOURCES ${source})
    endif()
endforeach()

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${HEADERS} ${SOURCES})

# =========================================================
# 2. 包含路径与输出设置 (Includes & Output)
# =========================================================
# 确保 ENGINE_ROOT_DIR 有值（防止单独打开此文件时路径为空）
if(NOT DEFINED ENGINE_ROOT_DIR)
    get_filename_component(ENGINE_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
endif()

set(LLVM_INCLUDE_DIRS 
    ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/LLVM/include 
    ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/mustache 
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# 输出目录统一设定到 Engine/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${ENGINE_ROOT_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${ENGINE_ROOT_DIR}/bin)

include_directories(${LLVM_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/parser)

# =========================================================
# 3. 目标构建 (Target Build)
# =========================================================
add_executable(${TARGET_NAME} ${HEADERS} ${SOURCES})

set_target_properties(${TARGET_NAME} PROPERTIES CXX_STANDARD 17)
set_target_properties(${TARGET_NAME} PROPERTIES FOLDER "Tools")

# 平台特定的链接设置
if (CMAKE_HOST_WIN32)
    set(LLVM_LIBRARY_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/LLVM/lib/x64)
    set(LLVM_SHARED_LIBRARY_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/LLVM/bin/x64)
    # 使用 -O2 而不是 /O2 以兼容 MinGW/GCC
    set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O2")
    target_link_libraries(${TARGET_NAME}  ${LLVM_LIBRARY_DIR}/libclang.lib)
elseif(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Linux")
    set(LLVM_LIBRARY_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/LLVM/lib/Linux)
    set(LLVM_SHARED_LIBRARY_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/LLVM/bin/Linux)
    set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O3")
    target_link_libraries(${TARGET_NAME}  ${LLVM_SHARED_LIBRARY_DIR}/libclang.so.12)
elseif(CMAKE_HOST_APPLE) 
    # macOS logic...
    set(LLVM_LIBRARY_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/LLVM/lib/macOS)
    set(LLVM_SHARED_LIBRARY_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/LLVM/bin/macOS)
    set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O3")
    target_link_libraries(${TARGET_NAME}  ${LLVM_SHARED_LIBRARY_DIR}/libclang.dylib)
endif()

add_definitions(-DTIXML_USE_STL)

# 构建后拷贝 libclang.dll 到 bin 目录
add_custom_command(TARGET ${TARGET_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${LLVM_SHARED_LIBRARY_DIR}"
        $<TARGET_FILE_DIR:${TARGET_NAME}>
)

# =========================================================
# 4. 预编译反射配置 (Reflection Setup)
# =========================================================
# 修正路径：工具生成的目录应该和 CMAKE_RUNTIME_OUTPUT_DIRECTORY 一致
set(PRECOMPILE_TOOLS_PATH "${ENGINE_ROOT_DIR}/bin") 

# 设置系统包含路径
if (CMAKE_HOST_WIN32)
    set(sys_include "*") 
elseif(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Linux" )
    set(sys_include "/usr/include/c++/9/") 
elseif(CMAKE_HOST_APPLE)
    find_program(XCRUN_EXECUTABLE xcrun)
    execute_process(
      COMMAND ${XCRUN_EXECUTABLE} --sdk macosx --show-sdk-platform-path
      OUTPUT_VARIABLE osx_sdk_platform_path_test
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(sys_include "${osx_sdk_platform_path_test}/../../Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1") 
endif()

# 聚合头文件路径 (通常生成在 Build 目录)
set(PARSER_INPUT "${CMAKE_BINARY_DIR}/parser_header.h")

# =========================================================
# 5. 执行反射生成 (Running Reflection)
# =========================================================
set(PRECOMPILE_TARGET "LIZERALPreCompile")

# 使用 Generator Expression 获取 exe 的绝对路径，避免手动拼接出错
set(PRECOMPILE_PARSER_CMD $<TARGET_FILE:${TARGET_NAME}>)

add_custom_target(${PRECOMPILE_TARGET} ALL
    COMMAND ${CMAKE_COMMAND} -E echo "*************************************************************"
    COMMAND ${CMAKE_COMMAND} -E echo "**** [Precompile] BEGIN "
    COMMAND ${CMAKE_COMMAND} -E echo "*************************************************************"
    
    # 运行解析器
    COMMAND ${PRECOMPILE_PARSER_CMD} 
        "${ENGINE_ROOT_DIR}/engine/source/runtime" 
        "${PARSER_INPUT}" 
        "${ENGINE_ROOT_DIR}/engine/source;${ENGINE_ROOT_DIR}/engine/3rdparty" 
        "${sys_include}" 
        "LIZERAL" 
        "0"
    
    COMMAND ${CMAKE_COMMAND} -E echo "+++ Precompile finished +++"
    
    # 打印提示，如果出错显示路径
    COMMENT "Running Reflection Tool: ${TARGET_NAME}"
)

# [关键修正] 添加依赖：确保先编译出 exe，再运行它
add_dependencies(${PRECOMPILE_TARGET} ${TARGET_NAME})

set_target_properties(${PRECOMPILE_TARGET} PROPERTIES FOLDER "Tools")