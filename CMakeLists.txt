cmake_minimum_required(VERSION 3.7)              # CMake version check
project(RESEARCH)                                # Create project "RESEARCH"

set(CMAKE_CXX_STANDARD 11)                       # Enable C++11 standard

set(SOURCE_FILES asylmssp.cpp)                   # Add asylmssp.cpp file of project root directory as source file

add_executable(MAIN main.cpp
        packing.cpp
        packing.h
        base.cpp
        base.h)                                  # Add executable target with source files listed in SOURCE_FILES variable
