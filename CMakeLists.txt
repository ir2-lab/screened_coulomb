cmake_minimum_required(VERSION 3.10.0)

project(screened_coulomb VERSION 0.1.0 LANGUAGES CXX)

project(screened_coulomb 
        VERSION 0.1.0
        DESCRIPTION "Header-only C++ library for Screened Coulomb scattering calculations"
        HOMEPAGE_URL https://gitlab.com/ir2-lab/
        LANGUAGES CXX)

##
## MAIN_PROJECT CHECK
## determine if it is built as a subproject (using add_subdirectory) or if it is the main project
##
set(MAIN_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(MAIN_PROJECT ON)
endif()

add_library(screened_coulomb INTERFACE)

target_include_directories(
    screened_coulomb
    INTERFACE
    include
)

if(MAIN_PROJECT)
  add_executable(test test/test.cpp)
  target_link_libraries(test screened_coulomb)
endif()

