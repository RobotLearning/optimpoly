cmake_minimum_required (VERSION 2.8)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wall -DNDEBUG -std=c++11 -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -DDEBUG -Wall -Wextra -pedantic -std=c++11 -pthread")

set(PROJECT_NAME table-tennis)
set(PROJECT_VERSION 1.2)
project(${PROJECT_NAME})

# INCLUDE ARMADILLO AS REQUIRED LIBRARY
#find_package(Armadillo REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake_files/cmake_findpkg")
#set (BOOST_ROOT $ENV{HOME}/install)
#set (ARMADILLO_LIBRARY $ENV{HOME}/install/lib)
#set (ARMADILLO_INCLUDE_DIR $ENV{HOME}/install/include)
#set (NLOPT_LIBRARIES $ENV{HOME}/install/lib)
#set (NLOPT_INCLUDE_DIRS $ENV{HOME}/install/include)
find_package (Boost COMPONENTS program_options unit_test_framework REQUIRED) 
# second one is only for testing so in the future test cmake could enforce only
find_package (Armadillo REQUIRED PATHS /home/robolab/software)
include_directories(${ARMADILLO_INCLUDE_DIRS})
find_package (NLOPT REQUIRED)
include_directories(${NLOPT_INCLUDE_DIRS})
find_library(ZMQPP NAMES zmqpp)
find_path(ZMQPP_INCLUDES NAMES zmqpp/zmqpp.hpp)

add_subdirectory(src)
add_subdirectory(test)
# our own library is tested also
#find_package(${PROJECT_NAME} REQUIRED)

# MAKE SURE TO ADD USER LIBRARY TO LD_LIBRARY_PATH!
