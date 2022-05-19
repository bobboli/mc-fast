add_executable(Profile "Main.cpp" )
add_executable(ProfileFromFile "ProfileFromFile.cpp")

target_link_libraries(Profile "MarchingCubes")
target_link_libraries(ProfileFromFile "MarchingCubes")
