# Bullet Physics Library CMakeLists.txt
# This file creates an interface library for Bullet Physics that provides
# convenient include paths for using bullet headers with #include "bullet/..."

cmake_minimum_required(VERSION 3.10)

# Enable policy CMP0079 to allow target_link_libraries with targets from other directories
cmake_policy(SET CMP0079 NEW)

# Set Bullet version variables if not already set
if(NOT DEFINED BULLET_VERSION)
    set(BULLET_VERSION "3.24")
endif()

if(NOT DEFINED BULLET_PHYSICS_SOURCE_DIR)
    set(BULLET_PHYSICS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
endif()

if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS OFF)
endif()

if(NOT DEFINED INSTALL_LIBS)
    set(INSTALL_LIBS OFF)
endif()

# Create an interface library for Bullet
add_library(bullet INTERFACE)

# Set the include directories for the bullet interface library
# This allows other targets to use #include "bullet/btBulletDynamicsCommon.h"
target_include_directories(bullet INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/BulletCollision
    ${CMAKE_CURRENT_SOURCE_DIR}/BulletDynamics
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearMath
)

# Also set include directories for the actual bullet libraries
# This ensures that when building bullet itself, the include paths are correct
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add subdirectories for the actual Bullet libraries
add_subdirectory(BulletCollision)
add_subdirectory(BulletDynamics)
add_subdirectory(LinearMath)

# Make the bullet interface library depend on the actual Bullet libraries
target_link_libraries(bullet INTERFACE
    BulletCollision
    BulletDynamics
    LinearMath
)

# Set up library dependencies between bullet components
target_link_libraries(BulletDynamics BulletCollision LinearMath)
target_link_libraries(BulletCollision LinearMath)

# Export the bullet target for use in other CMake projects
set_target_properties(bullet PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)

# Create an alias for backward compatibility
add_library(Bullet::Bullet ALIAS bullet)
