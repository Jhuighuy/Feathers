cmake_minimum_required(VERSION 3.18)
project(Feathers)

set(CMAKE_CXX_STANDARD 23)

include_directories(src)

file(GLOB_RECURSE SOURCES src/*.cc src/*.cxx src/*.hh src/*.hxx src/*.inl)

add_executable(Feathers ${SOURCES})

## ------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------ ##

target_include_directories(Feathers PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lib")

# Set path to custom modules.
list(APPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/cmake")

# Enable GLM.
find_package(GLM REQUIRED)

# Enable Intel TBB.
find_package(TBB REQUIRED tbb)
if(TBB_FOUND)
    add_compile_definitions(FEATHERS_HAS_TBB=1)
    target_link_libraries(Feathers ${TBB_IMPORTED_TARGETS})
else()
    # Fallback to OpenMP.
    find_package(OpenMP)
    if(OpenMP_FOUND)
        add_compile_definitions(FEATHERS_HAS_OPENMP=1)
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS
            "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

## ------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------ ##

set(CMAKE_CXX_FLAGS_DEBUG
    "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE
    "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -march=native -DNDEBUG=1")

## ------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------ ##
