cmake_minimum_required(VERSION 3.10)
project(oadcs_project)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math -funroll-loops")


# Add headers
include_directories(include)

# # Add astrolib
# include_directories(astrolib)
# include_directories(astrolib/include)

# Add Eigen
set(EIGEN3_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/external/eigen")
include_directories(${EIGEN3_INCLUDE_DIR})

# # Add raylib/raylib-cpp
# set(RAYLIB_DIR "${CMAKE_SOURCE_DIR}/external/raylib")
# include_directories("${RAYLIB_DIR}/src")
# set(RAYLIB_CPP_DIR "${CMAKE_SOURCE_DIR}/external/raylib-cpp")
# include_directories("${RAYLIB_CPP_DIR}/include")

# Add clay/clayman

# Automatically find all .cpp files
file(GLOB_RECURSE SOURCES 
    src/*.cpp 
    # astrolib/src/*.cpp
)

# Add executable
add_executable(oadcs_project ${SOURCES})

# Link raylib
add_subdirectory(${RAYLIB_DIR} EXCLUDE_FROM_ALL)
# target_link_libraries(oadcs_project raylib)
target_link_libraries(oadcs_project)