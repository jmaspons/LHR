cmake_minimum_required(VERSION 2.6)
project(r-package)

# cmake modules (compute path relative to this file)
get_filename_component(ROOT_SRC_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
set(CMAKE_MODULE_PATH "${ROOT_SRC_DIR}/cmake/modules/")

find_package(LibR REQUIRED)
include_directories(${LIBR_INCLUDE_DIRS})
add_library(betaBinomial SHARED betabinom.c betanbinom.c)

add_executable(r-package main.c)
message(STATUS "R: ${LIBR_LIBRARIES}\n")
target_link_libraries(r-package betaBinomial ${LIBR_LIBRARIES}) 
