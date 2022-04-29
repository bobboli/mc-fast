set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)  # todo: how to export only the necessary symbols?
add_library(MarchingCubes SHARED  "MarchingCubes.cpp" "mcTables.h" "MarchingCubes.h" "Geometry.h")
