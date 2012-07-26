cmake_minimum_required(VERSION 2.8)

project(LHR)

# cmake modules (compute path relative to this file)
get_filename_component(ROOT_SRC_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
set(CMAKE_MODULE_PATH "${ROOT_SRC_DIR}/src/cmake/modules/")
message(STATUS "Module path: ${CMAKE_MODULE_PATH}\n")

find_package(LibR REQUIRED)
include_directories(${LIBR_INCLUDE_DIRS})

add_subdirectory(src)
