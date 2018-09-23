cmake_minimum_required(VERSION 2.8)

file(GLOB cpp_files src/*.cpp)

add_library(maths_in_physic ${cpp_files})
target_include_directories(maths_in_physic PUBLIC include)
target_link_libraries(maths_in_physic PUBLIC armadillo)
