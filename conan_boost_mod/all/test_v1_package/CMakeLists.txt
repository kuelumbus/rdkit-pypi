cmake_minimum_required(VERSION 3.15)
project(test_v1_package LANGUAGES CXX)

enable_testing()

find_package(Boost REQUIRED CONFIG)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../test_package
                 ${CMAKE_CURRENT_BINARY_DIR}/test_package)
