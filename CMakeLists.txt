cmake_minimum_required(VERSION 3.19)
project(Feathers)

set(CMAKE_CXX_STANDARD 17)

include_directories(src)

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(Boost 1.74.0 REQUIRED COMPONENTS graph)

# find_package(GLM)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -DNDEBUG=1")

file(GLOB_RECURSE SOURCES src/*.cc)

add_executable(Feathers ${SOURCES})
