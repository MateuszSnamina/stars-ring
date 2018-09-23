cmake_minimum_required(VERSION 2.8)

file(GLOB cpp_files src/*.cpp)

add_library(stars_ring_basis ${cpp_files})
target_include_directories(stars_ring_basis PUBLIC include)
target_link_libraries(stars_ring_basis PUBLIC stars_ring_core maths_in_physic utility_kit armadillo)
