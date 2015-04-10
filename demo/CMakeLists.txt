cmake_minimum_required(VERSION 2.8)
add_executable(demoBiVariateGauss demoBiVariateGauss.cpp)
add_executable(demoPowSpecSamples2d demoPowSpecSamples2d.cpp)

file(COPY ${CMAKE_SOURCE_DIR}/demo/plotDemoBiVariateGauss.gp DESTINATION ${CMAKE_BINARY_DIR}/demo/)
file(COPY ${CMAKE_SOURCE_DIR}/demo/demoPowSpecSamples2d.gp DESTINATION ${CMAKE_BINARY_DIR}/demo/)
