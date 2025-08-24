cmake_minimum_required(VERSION 3.10.0)

project(screened_coulomb VERSION 0.1.0 LANGUAGES CXX)

project(screened_coulomb 
        VERSION 0.1.0
        DESCRIPTION "Header-only C++ library for Screened Coulomb scattering calculations"
        HOMEPAGE_URL https://gitlab.com/ir2-lab/
        LANGUAGES CXX)
        
# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(GNUInstallDirs)

add_library(screened_coulomb INTERFACE)

target_include_directories(
  screened_coulomb
  INTERFACE 
  ${CMAKE_INSTALL_INCLUDEDIR})

target_sources(screened_coulomb
  PUBLIC FILE_SET HEADERS
  BASE_DIRS ${PROJECT_SOURCE_DIR}/include
  FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/screened_coulomb.h
)

add_executable(test EXCLUDE_FROM_ALL test/test.cpp )
target_link_libraries(test screened_coulomb)

# install the headers that we specified with FILE_SET (header files for the user to include in their projects)
install(TARGETS screened_coulomb
  FILE_SET HEADERS
)