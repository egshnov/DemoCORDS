cmake_minimum_required(VERSION 3.26)
project(SoftFunctionalDependencies)

set(CMAKE_CXX_STANDARD 17)

add_executable(SoftFunctionalDependencies main.cpp
        candidate.h
        cords.hpp
        cords.cpp
        sample.hpp
        alliacies.hpp)
