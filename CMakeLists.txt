cmake_minimum_required(VERSION 3.5)
project(Rternatives)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -std=c++1y")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -std=c++1y -O3 -march=native")

include_directories(include)
set(SOURCE_FILES main.cpp include/CorTest.hxx include/utils.hxx include/Hypothesis.hxx)
add_executable(Rternatives ${SOURCE_FILES})