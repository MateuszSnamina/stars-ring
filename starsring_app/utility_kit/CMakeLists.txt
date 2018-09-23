cmake_minimum_required(VERSION 2.8)

file(GLOB cpp_files src/*.cpp)

add_library(utility_kit ${cpp_files})
target_include_directories(utility_kit PUBLIC include)
#target_link_libraries(utility_kit PUBLIC armadillo boost_regex)
