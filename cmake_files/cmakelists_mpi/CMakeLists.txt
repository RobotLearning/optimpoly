cmake_minimum_required (VERSION 2.8)

IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE RelWithDebInfo)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
set(WARNING_FLAGS "-Werror -Wall -Weffc++ -Wno-error=effc++ -Wextra -pedantic -pedantic-errors")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -DNDEBUG -std=c++11 -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG ${WARNING_FLAGS} -std=c++11 -pthread")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -Wall -std=c++11 -pthread")
set(BUILD_TEST FALSE CACHE BOOL "Build tests if ON")

set(PROJECT_NAME table-tennis)
set(PROJECT_VERSION 1.2)
project(${PROJECT_NAME})

# INCLUDE ARMADILLO AS REQUIRED LIBRARY
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake_files/cmake_findpkg")
find_package(Armadillo REQUIRED)
find_package (Boost COMPONENTS program_options unit_test_framework REQUIRED) 
include_directories(${ARMADILLO_INCLUDE_DIRS})
find_package (NLOPT REQUIRED)
include_directories(${NLOPT_INCLUDE_DIRS})
find_library(ZMQPP NAMES zmqpp)
find_path(ZMQPP_INCLUDES NAMES zmqpp/zmqpp.hpp)

add_subdirectory(src)
if (BUILD_TEST)
  add_subdirectory(test)
endif()



