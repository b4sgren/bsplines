cmake_minimum_required(VERSION 3.4)
project(bsplines CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED 1)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

option(BUILD_TESTS ON CACHE)

find_package(Eigen3 REQUIRED)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE
    include/bsplines/
)
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)

if (BUILD_TESTS)
    add_subdirectory(extern/googletest)
    add_subdirectory(tests)
endif()
