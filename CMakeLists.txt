cmake_minimum_required(VERSION 3.1)

project(BVH)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Change the build configurations available
set(CMAKE_CONFIGURATION_TYPES "Debug;Develop;Release" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_DEVELOP ${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_RELEASE ${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} CACHE STRING "" FORCE)

# Set tanks as the start up project
set_directory_properties(PROPERTIES VS_STARTUP_PROJECT BVHTest)

# ===========================================================
# TOOLS

set(PROJ_NAME BVHTest)
set(INCLUDE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/include/)
set(SOURCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/source/)

file(GLOB_RECURSE INCLUDE_FILES ${INCLUDE_ROOT}/*.*)
file(GLOB_RECURSE SOURCE_FILES ${SOURCE_ROOT}/*.*)
source_group(TREE ${INCLUDE_ROOT} FILES ${INCLUDE_FILES})
source_group(TREE ${SOURCE_ROOT} FILES ${SOURCE_FILES})
add_executable(${PROJ_NAME} ${INCLUDE_FILES} ${SOURCE_FILES})

target_include_directories(${PROJ_NAME} PUBLIC ${INCLUDE_ROOT})
target_include_directories(${PROJ_NAME} PUBLIC ${SOURCE_ROOT})
target_include_directories(${PROJ_NAME} PUBLIC ${SOURCE_ROOT}/imgui)

target_link_libraries(${PROJ_NAME} PUBLIC D3D11.lib) 
