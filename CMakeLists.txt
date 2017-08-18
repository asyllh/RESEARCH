cmake_minimum_required(VERSION 3.7)              # CMake version check
project(RESEARCH)                                # Create project "RESEARCH"

set(CMAKE_CXX_STANDARD 11)                       # Enable C++11 standard

set(SOURCE_FILES asylmssp.cpp)                   # Add asylmssp.cpp file of project root directory as source file

add_executable(AsylMSSP asylmssp.cpp)            # Add executable target with source files listed in SOURCE_FILES variable
add_executable(MSSPOutput msspoutput.cpp)
add_executable(Strips msspstrips.cpp)
add_executable(MAIN main.cpp mbahra.cpp mbahra.h packing.cpp packing.h)