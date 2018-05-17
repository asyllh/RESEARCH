cmake_minimum_required(VERSION 3.8)
project(COMBINED)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp packing.cpp packing.h base.cpp base.h)
add_executable(COMBINED ${SOURCE_FILES})