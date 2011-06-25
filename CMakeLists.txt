project(cmpt885-hw2)

cmake_minimum_required(VERSION 2.8.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp")

add_definitions(-Wall -Wextra)

add_subdirectory(genetic)
