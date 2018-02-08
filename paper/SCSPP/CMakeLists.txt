cmake_minimum_required(VERSION 3.8)
project(ffdexact)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        base.cpp
        base.h
        main.cpp
        packing.cpp
        packing.h)

add_executable(ffdexact ${SOURCE_FILES})