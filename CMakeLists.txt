cmake_minimum_required(VERSION 3.5)
project(Rternatives)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

include_directories(include)
set(SOURCE_FILES main.cpp include/CorTest.hxx include/utils.hxx)
add_executable(Rternatives ${SOURCE_FILES})