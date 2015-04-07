cmake_minimum_required(VERSION 2.8)

# Guard against in-source builds (got this from Eigen)
IF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
ENDIF()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

FIND_PACKAGE(Eigen3 REQUIRED)
find_package(nanoflann REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${NANOFLANN_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -Wextra -Wfatal-errors")

ENABLE_TESTING()

add_subdirectory(test)
add_subdirectory(bench)
add_subdirectory(data)
add_subdirectory(demo)
