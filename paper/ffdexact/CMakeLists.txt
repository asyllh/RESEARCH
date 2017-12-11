cmake_minimum_required(VERSION 3.8)
project(ffdexact)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        baseexact.cpp
        baseexact.h
        mainexact.cpp
        packingexact.cpp
        packingexact.h)

add_executable(ffdexact ${SOURCE_FILES})