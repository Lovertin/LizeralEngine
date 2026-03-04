cmake_minimum_required(VERSION 3.16)
project(LizeralEngine VERSION 0.0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#导入qt6本地文件地址
set(CMAKE_PREFIX_PATH "C:/Qt/6.10.0/mingw_64" ${CMAKE_PREFIX_PATH})

find_package(Qt6 COMPONENTS Core Gui Widgets OpenGLWidgets REQUIRED)

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

# A. 核心数学库 (Math)
file(GLOB_RECURSE SRC_MATH "source/runtime/core/math/*.cpp")

# B. ECS 系统
file(GLOB_RECURSE SRC_ECS "source/runtime/function/ecs/*.cpp")

# C. 物理系统
file(GLOB_RECURSE SRC_PHYSICS "source/runtime/function/physics/*.cpp")

# D. 必要的 Log 或其他 Core (根据你的项目结构调整，只选纯 cpp)
file(GLOB_RECURSE SRC_LOG "source/runtime/core/log/*.cpp")

# 新增 G: 抓取 Meta 系统 (反射/序列化底层核心实现) 的源文件
file(GLOB_RECURSE SRC_META "source/runtime/core/meta/*.cpp")

# 新增 H: 抓取由 LizeralParser 自动生成的所有 .cpp 文件
file(GLOB_RECURSE SRC_GENERATED "source/_generated/*.cpp")

# 1. 确保 Input 源文件被包含
file(GLOB_RECURSE SRC_INPUT "source/runtime/function/input/*.cpp") 

# 2. 确保 Render (Camera) 源文件被包含
file(GLOB_RECURSE SRC_RENDER "source/runtime/function/render/*.cpp")

# 3. res_type 系统 (用于资源管理)
file(GLOB_RECURSE SRC_TYPE "source/runtime/function/res_type/*.cpp")

# 4. Resource 系统 (用于资源管理)
file(GLOB_RECURSE SRC_RESOURCE "source/runtime/resource/*.cpp")

# 5. Vendor包
file(GLOB_RECURSE SRC_VENFOR "source/runtime/core/vendor/stb_image/*.cpp")

# 6. 包含我的editor
file(GLOB_RECURSE SRC_EDITOR CONFIGURE_DEPENDS "source/editor/*.cpp" "source/editor/*.h")



# E. 聚合物理沙盒源文件
set(PHYSICS_SANDBOX_SOURCES 
    ${SRC_MATH}
    ${SRC_ECS}
    ${SRC_PHYSICS}
    ${SRC_LOG}
    ${SRC_INPUT}
    ${SRC_RENDER}
    ${SRC_TYPE}
    ${SRC_RESOURCE}
    ${SRC_VENFOR}
    "source/runtime/function/sandbox/PhysicsBaseTest/sandbox.cpp"
)

set(EDITOR_SOURCES 
    ${SRC_MATH}
    ${SRC_ECS}
    ${SRC_PHYSICS}
    ${SRC_LOG}
    ${SRC_INPUT}  
    ${SRC_RENDER} 
    ${SRC_TYPE}
    ${SRC_RESOURCE}
    ${SRC_VENFOR}
    ${SRC_META}      
    ${SRC_GENERATED} 
    "3rdparty/json11/json11.cpp"
    ${SRC_EDITOR}
)



# =============================================================
# ★ 本地引入工业级模型加载库：Assimp
# =============================================================
# 关闭 Assimp 那些我们不需要的附属工具和测试，极大加快编译速度
set(ASSIMP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_ASSIMP_TOOLS OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_ZLIB ON CACHE BOOL "" FORCE) # 确保开启 zlib 支持
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE) # 静态链接，免去处理 DLL 的烦恼

# 直接将本地的 assimp 文件夹加入编译树
add_subdirectory(source/runtime/core/vendor/assimp)

# 包含 Assimp 的全局头文件路径，这样你在代码里才能 #include <assimp/Importer.hpp>
include_directories(source/runtime/core/vendor/assimp/include)

# 针对 assimp 目标禁用未初始化变量警告
if(MINGW OR (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "13.0"))
    # 检查 assimp 的实际目标名称
    if(TARGET assimp)
        target_compile_options(assimp PRIVATE 
            -Wno-error=maybe-uninitialized
            -Wno-maybe-uninitialized
        )
        message(STATUS "Applied GCC 13+ compatibility fixes to assimp")
    elseif(TARGET assimp::assimp)
        # 如果是别名目标，不能直接使用 target_compile_options
        # 改为使用 set_target_properties 为源文件设置属性
        get_target_property(ASSIMP_SOURCES assimp::assimp SOURCES)
        if(ASSIMP_SOURCES)
            set_source_files_properties(${ASSIMP_SOURCES} PROPERTIES
                COMPILE_FLAGS "-Wno-error=maybe-uninitialized -Wno-maybe-uninitialized"
            )
            message(STATUS "Applied GCC 13+ compatibility fixes to assimp sources")
        endif()
    endif()
endif()


set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

#构建qt6平台
add_executable(LizeralEditor ${EDITOR_SOURCES})
target_include_directories(LizeralEditor PRIVATE
    source                      
    source/runtime              
    3rdparty/glfw/include       
    3rdparty/json11             
)

target_link_libraries(LizeralEditor PUBLIC
    glad        
    glfw
    bullet
    assimp      
    # ★ 新增：链接 Qt 模块
    Qt6::Core 
    Qt6::Gui 
    Qt6::Widgets 
    Qt6::OpenGLWidgets
)

if(WIN32)
    target_link_libraries(LizeralEditor PUBLIC opengl32)
endif()


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
)

# 链接库
target_link_libraries(LizeralPhysicsSandbox PUBLIC
    glad        # 【关键】链接 glad 库，自动修复 glad.h 和 khrplatform.h 问题
    glfw
    bullet
    assimp      # 【关键】链接本地编译好的 Assimp
)

if(WIN32)
    target_link_libraries(LizeralPhysicsSandbox PUBLIC opengl32)
endif()

# =============================================================
# 5. 构建渲染测试程序
# =============================================================
# add_executable(LizeralRenderingTest ${RENDERING_TEST_SOURCES})

# # 设置包含路径 (Header Search Paths)
# target_include_directories(LizeralRenderingTest PRIVATE
#     source                      # 允许 #include "runtime/..."
#     source/runtime              # 允许 #include "core/..."
#     3rdparty/glfw/include       # GLFW
#     3rdparty/json11             # json11.hpp
# )

# # 链接库
# target_link_libraries(LizeralRenderingTest PUBLIC
#     glad        
#     glfw
#     assimp      # 【关键】为渲染测试也链接上 Assimp，方便后续测试模型
# )

# if(WIN32)
#     target_link_libraries(LizeralRenderingTest PUBLIC opengl32)
# endif()