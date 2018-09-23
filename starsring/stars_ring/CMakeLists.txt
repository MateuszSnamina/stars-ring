cmake_minimum_required(VERSION 2.8)

file(GLOB cpp_files src/*.cpp)

add_executable(stars_ring ${cpp_files})
target_include_directories(stars_ring PUBLIC include)
target_link_libraries(stars_ring PUBLIC utility_kit maths_in_physic maths_in_physic stars_ring_basis stars_ring_numerical stars_ring_analytical armadillo boost_program_options)
