cmake_minimum_required(VERSION 3.0)

project(xRS2)

include_directories(${CMAKE_SOURCE_DIR}/include)

SET( CMAKE_CXX_FLAGS "-std=c++11 -Wall")

SET( CMAKE_BUILD_TYPE Debug)

add_executable(matching_main src/matching_main.cc src/matching_generator.cc src/stripe.cc)