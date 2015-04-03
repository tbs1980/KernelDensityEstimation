cmake_minimum_required(VERSION 2.8)

project(testKDE)

#add_executable(testKDE testKDE.cpp)
add_executable(testGaussian Gaussian.cpp)
add_executable(testDiagonalBandwidthMatrix DiagonalBandwidthMatrix.cpp)

add_test(kernel_Gaussian testGaussian)
add_test(bandwidth_DiagonalBandwidthMatrix DiagonalBandwidthMatrix)
