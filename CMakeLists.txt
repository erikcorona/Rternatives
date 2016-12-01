cmake_minimum_required(VERSION 3.5)
project(Rternatives)

find_package(Boost COMPONENTS system filesystem thread REQUIRED)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall -std=c++1y")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -std=c++1y -O3 -march=native")

include_directories(include)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

if (APPLE OR WIN32)
    set(Boost_USE_STATIC_LIBS TRUE)
endif()

#set(CMAKE_CXX_COMPILER /usr/bin/clang++)
set(SOURCE_FILES main.cpp include/CorTest.hxx include/utils.hxx include/Hypothesis.hxx include/CommonOperations.hxx)
add_executable(Rternatives ${SOURCE_FILES})

add_subdirectory(tests)

#...
target_link_libraries(Rternatives
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_LIBRARIES}
        )