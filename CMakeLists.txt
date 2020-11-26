cmake_minimum_required(VERSION 3.12)

project(JacobiSVD)
add_executable(JacobiSVD main.cpp src/utils/util.cpp src/utils/matrix.cpp 
src/utils/block.cpp src/svd/two-sided/svd_blocked.cpp src/svd/two-sided/svd_subprocedure.cpp
src/svd/two-sided/nsvd.cpp)

# Includes
#include(CTest)

# Set C++ standard and compiler flags.
set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-write-strings -DNDEBUG -O3 -mavx2 -mfma -ftree-vectorize")

# If cmake is not passed the CMAKE_BUILD_TYPE=Debug option, default to Release.
# Note that if the CMAKE_BUILD_TYPE option is passed when running cmake, its value is cached
# for subsequent runs. To clear it delete the cmake cache or set it to another build type.
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(Yellow      "${Esc}[33m")
endif()

# Print the build type in yellow color to make it visible.
message("${Yellow}Build type set to ${CMAKE_BUILD_TYPE}${ColourReset}")

# Project subdirectories
add_subdirectory(src)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()