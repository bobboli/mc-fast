# CMakeList.txt : CMake project for mc-opt, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)
project(mc-opt)


# Add source to this project's executable.

# TODO: Add tests and install targets if needed.

add_subdirectory("MarchingCubes")
add_subdirectory("Profile")

