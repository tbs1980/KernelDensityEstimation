cmake_minimum_required(VERSION 2.8)
add_executable(demoBiVariateGauss demoBiVariateGauss.cpp)

file(COPY ${CMAKE_SOURCE_DIR}/demo/plotDemoBiVariateGauss.gp DESTINATION ${CMAKE_BINARY_DIR}/demo/)
