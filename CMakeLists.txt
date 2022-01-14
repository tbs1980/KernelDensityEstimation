cmake_minimum_required(VERSION 3.1)

project(
  KernelDensityEstimation
  VERSION "1.0.0"
  DESCRIPTION ""
  LANGUAGES C CXX)

include(FetchContent)

option(ENABLE_BUILD_TEST "Build the tests of the project" OFF)
option(ENABLE_BUILD_BENCHMARKS "Build the tests of the project" OFF)
option(ENABLE_BUILD_DEMO "Build the tests of the project" OFF)

# Guard against in-source builds (got this from Eigen)
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. "
  )
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

FetchContent_Declare(
  nanoflann
  GIT_REPOSITORY "https://github.com/jlblancoc/nanoflann"
  GIT_TAG "b3e81831b847a77473e0da2ad7ee266b4f4e0d9a")

# Get nanoflan lib
find_package(nanoflann)
if(NOT nanoflann_FOUND)
  message("Fetching ynanoflann lib ...")
  FetchContent_MakeAvailable(nanoflann)
endif()

find_package(Eigen3 REQUIRED)

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -Wextra -Wfatal-errors")

if(ENABLE_BUILD_TEST)
  enable_testing()
  add_subdirectory(test)
endif(ENABLE_BUILD_TEST)

if(ENABLE_BUILD_BENCHMARKS)
  add_subdirectory(bench)
endif(ENABLE_BUILD_BENCHMARKS)

if(ENABLE_BUILD_DEMO)
  add_subdirectory(data)
  add_subdirectory(demo)
endif(ENABLE_BUILD_DEMO)

add_library(KernelDensityEstimation INTERFACE)
target_include_directories(KernelDensityEstimation INTERFACE src/)
target_link_libraries(KernelDensityEstimation INTERFACE Eigen3 nanoflann)
