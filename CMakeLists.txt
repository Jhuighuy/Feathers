cmake_minimum_required(VERSION 3.19)
project(Feathers)

set(CMAKE_CXX_STANDARD 17)

include_directories(src)

## ------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------ ##

# Enable OpenMP for release builds.
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS_RELEASE
        "${CMAKE_C_FLAGS_RELEASE} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS_RELEASE
        "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE
        "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# find_package(GLM)

## ------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------ ##

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -DNDEBUG=1")

## ------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------ ##

file(GLOB_RECURSE SOURCES src/*.cc)

add_executable(Feathers ${SOURCES})
