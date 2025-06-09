cmake_minimum_required(VERSION 3.31)
project(SCCG)

set(CMAKE_CXX_STANDARD 20)

add_executable(SCCG main.cpp src/Compressor.cpp src/KmerIndex.cpp src/MatchFinder.cpp src/Position.h src/Utils.h)
