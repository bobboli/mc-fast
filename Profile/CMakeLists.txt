add_executable(Profile "ProfileSphere.cpp" "SimpleTimer.cpp" "SimplexNoise.cpp")
add_executable(ProfileFromFile "ProfileFromFile.cpp" "SimpleTimer.cpp")

target_link_libraries(Profile "MarchingCubes")
target_link_libraries(ProfileFromFile "MarchingCubes")
