cmake_minimum_required(VERSION 3.31)
project(SCCG)

set(CMAKE_CXX_STANDARD 20)

add_executable(SCCG main.cpp
        src/KmerIndex.cpp
        src/KmerIndex.h
        src/Utils.h
        src/MatchFinder.cpp
        src/MatchFinder.h
        src/Compressor.cpp
        src/Compressor.h
)
