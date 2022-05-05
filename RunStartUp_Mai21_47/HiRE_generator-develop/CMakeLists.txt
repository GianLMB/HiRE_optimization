cmake_minimum_required(VERSION 2.6)
project(opep_cpp)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -O2")

set(mainlib_sources functions.cpp fa2cg.cpp)
add_library(mainlib ${mainlib_sources})

add_subdirectory(entities)
add_subdirectory(geometry)
add_subdirectory(inout)
add_subdirectory(log)
add_subdirectory(mathematics)
add_subdirectory(param)

set(program HiRE_parm)
set(HiRElibs entities geometry inout log mathematics param)

add_executable(${program} main.cpp)
add_executable(fa2cg main_fa2cg.cpp)
target_link_libraries(${program} mainlib ${HiRElibs})
target_link_libraries(fa2cg mainlib ${HiRElibs})

