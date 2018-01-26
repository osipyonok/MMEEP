cmake_minimum_required(VERSION 3.9)
project(MEEP_1)

set(CMAKE_CXX_STANDARD 17)

add_executable(MEEP_1 main.cpp minpack.cpp minpack.h gnuplot_i.hpp)