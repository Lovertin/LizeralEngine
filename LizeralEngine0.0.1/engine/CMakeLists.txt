cmake_minimum_required(VERSION 3.16)
project(LizeralEngine VERSION 0.0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 二进制输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/bin)

# =============================================================
# 1. 解决 GLAD 问题：将其封装为静态库 (Target)
# =============================================================
# 这是一个现代 CMake 技巧。我们创建一个名为 'glad' 的库目标。
# 任何链接到 'glad' 的项目，都会自动继承它的 include 路径。
add_library(glad STATIC "3rdparty/glad/src/glad.c")

# 【关键】把 include 路径“绑定”到 glad 库上
# 需要包含 glad/include 和 glad/include/KHR 目录
target_include_directories(glad PUBLIC 
    "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/glad/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/glad/include/KHR"
)

# =============================================================
# 2. 第三方库配置
# =============================================================
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)

add_subdirectory(3rdparty/glfw)
add_subdirectory(3rdparty/bullet)

# =============================================================
# 3. 解决 static_assert 问题：白名单选取源码
# =============================================================
# 我们不再 GLOB_RECURSE 整个 source 目录，而是只抓取物理引擎真正需要的。
# 这样可以绝对避开 'meta/serializer' 和 'tool/parser' 里的报错代码。

# A. 核心数学库 (Math)
file(GLOB_RECURSE SRC_MATH "source/runtime/core/math/*.cpp")

# B. ECS 系统
file(GLOB_RECURSE SRC_ECS "source/runtime/function/ecs/*.cpp")

# C. 物理系统
file(GLOB_RECURSE SRC_PHYSICS "source/runtime/function/physics/*.cpp")

# D. 必要的 Log 或其他 Core (根据你的项目结构调整，只选纯 cpp)
file(GLOB_RECURSE SRC_LOG "source/runtime/core/log/*.cpp")

# 1. 确保 Input 源文件被包含
# 假设你的目录结构是 source/runtime/function/input/input.cpp
file(GLOB_RECURSE SRC_INPUT "source/runtime/function/input/*.cpp") 

# 2. 确保 Render (Camera) 源文件被包含
# 假设你的目录结构是 source/runtime/function/render/CameraSystem.cpp 等
file(GLOB_RECURSE SRC_RENDER "source/runtime/function/render/*.cpp")

# 3. Pipeline 系统 (用于渲染测试)
file(GLOB_RECURSE SRC_PIPELINE "source/runtime/function/pipeline/*.cpp")

# 4. Resource 系统 (用于资源管理)
file(GLOB_RECURSE SRC_RESOURCE "source/runtime/resource/*.cpp")

file(GLOB_RECURSE SRC_MATERIAL "source/runtime/function/pipeline/Material/*.h")


# E. 聚合物理沙盒源文件
set(PHYSICS_SANDBOX_SOURCES 
    ${SRC_MATH}
    ${SRC_ECS}
    ${SRC_PHYSICS}
    ${SRC_LOG}
    ${SRC_INPUT}   # <--- 必须加这个
    ${SRC_RENDER}  # <--- 必须加这个
    ${SRC_PIPELINE}
    ${SRC_RESOURCE}
    ${SRC_MATERIAL}
    "source/runtime/function/sandbox/PhysicsBaseTest/sandbox.cpp"
)

# F. 聚合渲染测试源文件
set(RENDERING_TEST_SOURCES 
    ${SRC_MATH}
    ${SRC_ECS}
    ${SRC_LOG}
    ${SRC_INPUT}
    ${SRC_RENDER}
    ${SRC_PIPELINE}
    ${SRC_RESOURCE}
    "source/runtime/function/sandbox/renderingWorld.cpp"
)

# =============================================================
# 4. 构建物理沙盒
# =============================================================
add_executable(LizeralPhysicsSandbox ${PHYSICS_SANDBOX_SOURCES})

# 设置包含路径 (Header Search Paths)
target_include_directories(LizeralPhysicsSandbox PRIVATE
    source                      # 允许 #include "runtime/..."
    source/runtime              # 允许 #include "core/..."
    3rdparty/glfw/include       # GLFW
    3rdparty/json11             # json11.hpp
    # 注意：GLAD 的路径不需要在这里写，因为下面链接 glad 库时会自动带过来
)

# 链接库
target_link_libraries(LizeralPhysicsSandbox PUBLIC
    glad        # 【关键】链接 glad 库，自动修复 glad.h 和 khrplatform.h 问题
    glfw
    bullet
)

if(WIN32)
    target_link_libraries(LizeralPhysicsSandbox PUBLIC opengl32)
endif()

# =============================================================
# 5. 构建渲染测试程序
# =============================================================
add_executable(LizeralRenderingTest ${RENDERING_TEST_SOURCES})

# 设置包含路径 (Header Search Paths)
target_include_directories(LizeralRenderingTest PRIVATE
    source                      # 允许 #include "runtime/..."
    source/runtime              # 允许 #include "core/..."
    3rdparty/glfw/include       # GLFW
    3rdparty/json11             # json11.hpp
    # 注意：GLAD 的路径不需要在这里写，因为下面链接 glad 库时会自动带过来
)

# 链接库
target_link_libraries(LizeralRenderingTest PUBLIC
    glad        # 【关键】链接 glad 库，自动修复 glad.h 和 khrplatform.h 问题
    glfw
)

if(WIN32)
    target_link_libraries(LizeralRenderingTest PUBLIC opengl32)
endif()
