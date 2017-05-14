project(UnscentedKF)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++0x)

set(sources
    src/UKF.cpp
    src/main.cpp
    src/Tools.cpp)

add_executable(UnscentedKF ${sources})
