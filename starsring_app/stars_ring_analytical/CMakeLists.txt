cmake_minimum_required(VERSION 2.8)

file(GLOB cpp_files src/*.cpp)

add_library(stars_ring_analytical ${cpp_files})
target_include_directories(stars_ring_analytical PUBLIC include)
target_link_libraries(stars_ring_analytical PUBLIC stars_ring_core armadillo)
