# CMakeList.txt : CMake project for mc-opt, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)
project(mc-opt)


# Add source to this project's executable.

# TODO: Add tests and install targets if needed.

add_subdirectory("MarchingCubes")
add_subdirectory("Profile")

message(STATUS "${CMAKE_C_COMPILER}")   
message(STATUS "CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")   
message(STATUS "CMAKE_C_FLAGS_DEBUG = ${CMAKE_C_FLAGS_DEBUG}")   
message(STATUS "CMAKE_C_FLAGS_RELEASE = ${CMAKE_C_FLAGS_RELEASE}")   
message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")   
message(STATUS "CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")   
message(STATUS "CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}")  