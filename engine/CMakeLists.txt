cmake_minimum_required(VERSION 3.16)
project(LizeralEngine VERSION 0.0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_PREFIX_PATH "C:/Qt/6.10.0/mingw_64" ${CMAKE_PREFIX_PATH})

find_package(Qt6 COMPONENTS Core Gui Widgets OpenGLWidgets REQUIRED)

find_package(Vulkan REQUIRED)
if(Vulkan_FOUND)
    message(STATUS "====== Vulkan SDK Found! ======")
    message(STATUS "Vulkan Include Dir: ${Vulkan_INCLUDE_DIRS}")
    message(STATUS "Vulkan Library: ${Vulkan_LIBRARIES}")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/bin)

add_library(glad STATIC "3rdparty/glad/src/glad.c")

target_include_directories(glad PUBLIC 
    "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/glad/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/glad/include/KHR"
)

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)

add_subdirectory(3rdparty/glfw)
add_subdirectory(3rdparty/bullet)

# A. Math
file(GLOB_RECURSE SRC_MATH "source/runtime/core/math/*.cpp")

# B. ECS System
file(GLOB_RECURSE SRC_ECS "source/runtime/function/ecs/*.cpp")

# C. Physics System
file(GLOB_RECURSE SRC_PHYSICS "source/runtime/function/physics/*.cpp")

file(GLOB_RECURSE SRC_LOG "source/runtime/core/log/*.cpp")

# D. Reflection
file(GLOB_RECURSE SRC_META "source/runtime/core/meta/*.cpp")

# E. Generated
file(GLOB_RECURSE SRC_GENERATED "source/_generated/*.cpp")

# F. Input 
file(GLOB_RECURSE SRC_INPUT "source/runtime/function/input/*.cpp") 

# G. Rendering
file(GLOB_RECURSE SRC_RENDER "source/runtime/function/render/*.cpp")

# H. res_type 
file(GLOB_RECURSE VULKAN_SRC_TYPE "source/runtime/function/Vulkan_res_type/*.cpp")

file(GLOB_RECURSE SRC_TYPE "source/runtime/function/res_type/*.cpp")

# I. Resource 
file(GLOB_RECURSE SRC_RESOURCE CONFIGURE_DEPENDS "source/runtime/resource/*.cpp")

# J. Vendor package
file(GLOB_RECURSE SRC_VENFOR "source/runtime/core/vendor/stb_image/*.cpp")

# K. Editor
file(GLOB_RECURSE SRC_EDITOR CONFIGURE_DEPENDS "source/editor/*.cpp" "source/editor/*.h")

# L. meshopitimizer
file(GLOB MESHOPTIMIZER_SOURCES "3rdparty/meshopitimizer/src/*.cpp")


add_library(meshoptimizer STATIC ${MESHOPTIMIZER_SOURCES})


target_include_directories(meshoptimizer PUBLIC "3rdparty/meshopitimizer/src")


set(PHYSICS_SANDBOX_SOURCES 
    ${SRC_MATH}
    ${SRC_ECS}
    ${SRC_PHYSICS}
    ${SRC_LOG}
    ${SRC_INPUT}  
    ${SRC_RENDER} 
    ${SRC_TYPE}
    ${VULKAN_SRC_TYPE}
    ${SRC_RESOURCE}
    ${SRC_VENFOR}
    ${SRC_META}      
    ${SRC_GENERATED} 
    "3rdparty/json11/json11.cpp"
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
    ${VULKAN_SRC_TYPE}
    ${SRC_RESOURCE}
    ${SRC_VENFOR}
    ${SRC_META}      
    ${SRC_GENERATED} 
    "3rdparty/json11/json11.cpp"
    ${SRC_EDITOR}
)


set(ASSIMP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_ASSIMP_TOOLS OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_ZLIB ON CACHE BOOL "" FORCE)
set(ASSIMP_WARNINGS_AS_ERRORS OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_ALL_IMPORTERS_BY_DEFAULT OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_ALL_EXPORTERS_BY_DEFAULT OFF CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_GLTF_IMPORTER ON CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_OBJ_IMPORTER ON CACHE BOOL "" FORCE)
set(ASSIMP_BUILD_FBX_IMPORTER ON CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE) 

add_subdirectory(source/runtime/core/vendor/assimp)

include_directories(source/runtime/core/vendor/assimp/include)

if(MINGW OR (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "13.0"))
    if(TARGET assimp)
        target_compile_options(assimp PRIVATE 
            -Wno-error=maybe-uninitialized
            -Wno-maybe-uninitialized
            -Wno-error=unknown-pragmas
            -Wno-unknown-pragmas
        )
        message(STATUS "Applied GCC 13+ compatibility fixes to assimp")
    elseif(TARGET assimp::assimp)
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

add_compile_definitions(
    LIZERAL_ENGINE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}"
    LIZERAL_PROJECT_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/.."
    LIZERAL_SHADER_DIR="${CMAKE_CURRENT_SOURCE_DIR}/source/shader/"
)

# build qt6
add_executable(LizeralEditor ${EDITOR_SOURCES})
target_include_directories(LizeralEditor PRIVATE
    source                      
    source/runtime              
    3rdparty/glfw/include       
    3rdparty/json11   
    3rdparty/vma/include
    ${Vulkan_INCLUDE_DIRS}      
)

target_link_libraries(LizeralEditor PUBLIC
    glad        
    glfw
    bullet
    assimp      
    Qt6::Core 
    Qt6::Gui 
    Qt6::Widgets 
    Qt6::OpenGLWidgets
    Vulkan::Vulkan
    meshoptimizer
)

if(WIN32)
    target_link_libraries(LizeralEditor PUBLIC opengl32)
endif()

add_executable(LizeralPhysicsSandbox ${PHYSICS_SANDBOX_SOURCES})

target_include_directories(LizeralPhysicsSandbox PRIVATE
    source                     
    source/runtime             
    3rdparty/glfw/include      
    3rdparty/json11           
    3rdparty/vma/include
    ${Vulkan_INCLUDE_DIRS}
)

target_link_libraries(LizeralPhysicsSandbox PUBLIC
    glad        
    glfw
    bullet
    assimp
    Qt6::Core 
    Qt6::Gui 
    Qt6::Widgets 
    Qt6::OpenGLWidgets
    Vulkan::Vulkan
    meshoptimizer
)

if(WIN32)
    target_link_libraries(LizeralPhysicsSandbox PUBLIC opengl32)
endif()


# add_executable(LizeralRenderingTest ${RENDERING_TEST_SOURCES})

# # 设置包含路径 (Header Search Paths)
# target_include_directories(LizeralRenderingTest PRIVATE
#     source                     
#     source/runtime             
#     3rdparty/glfw/include      
#     3rdparty/json11         
# )

# target_link_libraries(LizeralRenderingTest PUBLIC
#     glad        
#     glfw
#     assimp     
# )

# if(WIN32)
#     target_link_libraries(LizeralRenderingTest PUBLIC opengl32)
# endif()
