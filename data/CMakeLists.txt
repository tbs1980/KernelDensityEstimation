cmake_minimum_required(VERSION 2.8)
add_executable(genMultiVariateGaussian genMultiVariateGaussian.cpp)

#copy test data to bin
file(COPY ${CMAKE_SOURCE_DIR}/data/powSpecSamples2D.dat DESTINATION ${CMAKE_BINARY_DIR}/data/)
file(COPY ${CMAKE_SOURCE_DIR}/data/bi-variate-Gauss.dat DESTINATION ${CMAKE_BINARY_DIR}/data/)
