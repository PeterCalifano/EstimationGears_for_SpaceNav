cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(estimation_gears_consumer VERSION 0.1.0 LANGUAGES CXX)

find_package(EstimationGears_for_SpaceNav REQUIRED CONFIG)

# Add executable target
add_executable(estimation_gears_consumer example_project.cpp)

target_link_libraries(
    estimation_gears_consumer
    PRIVATE EstimationGears_for_SpaceNav::EstimationGears_for_SpaceNav)
