cmake_minimum_required (VERSION 2.8)
project(${PROJECT_NAME})

find_library(ZMQ zmq)
# Using BOOST UNIT TEST FRAMEWORK for testing
#ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK) 

# CREATE SHARED LIBRARY
set(TEST_EXEC unit_tests)
set(TEST_SRC
    test_kalman.cpp
    test_kinematics.cpp
    test_optim.cpp
    test_table_tennis.cpp
    test_serve.cpp
    test_interface.cpp
)
add_executable (${TEST_EXEC} ${TEST_SRC})

# SET PROJECT VERSION
set_target_properties(${TEST_EXEC} PROPERTIES 
    VERSION ${PROJECT_VERSION})

# INCLUDE HEADERS (top folder, as opposed to CMAKE_CURRENT_SOURCE_DIR)
target_include_directories (${TEST_EXEC} PRIVATE
    ${ZMQPP_INCLUDES}
    ${CMAKE_SOURCE_DIR}/include #for json.hpp
    ${CMAKE_SOURCE_DIR}/include/serve
    ${CMAKE_SOURCE_DIR}/include/optim
    ${CMAKE_SOURCE_DIR}/include/player)

# BOOST UNIT TEST should be included statically whenever
# we use only init_unit_test_suite syntax!
#message("suffix: ${CMAKE_FIND_LIBRARY_SUFFIXES}")
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#message("suffix: ${CMAKE_FIND_LIBRARY_SUFFIXES}")
#find_library(UNIT_TEST NAMES libboost_unit_test_framework.a)
#MESSAGE(STATUS "BOOST UNIT TEST LOCATION: " ${UNIT_TEST})

add_library(boost_unit_test_framework STATIC IMPORTED)
set_target_properties(boost_unit_test_framework PROPERTIES IMPORTED_LOCATION $ENV{HOME}/install/lib/libboost_unit_test_framework.a)
#set_target_properties(boost_unit_test_framework PROPERTIES IMPORTED_LOCATION /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.a)
#set_target_properties(boost_unit_test_framework PROPERTIES IMPORTED_LOCATION /home/robolab/software/lib/libboost_unit_test_framework.a)
#add_library(boost_program_options SHARED IMPORTED)
#set_target_properties(boost_program_options PROPERTIES IMPORTED_LOCATION /home/robolab/software/lib/libboost_program_options.so)
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")

# INCLUDE "OUR" SHARED LIBRARY
target_link_libraries(${TEST_EXEC}
    ${ZMQPP}
    ${ZMQ}
    ${PROJECT_NAME}
    armadillo
    boost_program_options
    boost_unit_test_framework
    nlopt)

# INSTALL FOLDER
#find_library(${SHARED_OBJ} ${PROJECT_NAME} HINTS ${CMAKE_SOURCE_DIR}/lib)
#target_link_libraries(${TEST_EXEC} ${SHARED_OBJ})
install(TARGETS ${TEST_EXEC}
    DESTINATION ${CMAKE_SOURCE_DIR})

# OPTIMIZATION OPTIONS
#target_compile_options(${TEST_EXEC} PRIVATE
#                -O3 -Wall -DNDEBUG -std=c++11 -pthread)

# DEBUG OPTIONS
#target_compile_options(example PUBLIC 
#        -g -DDEBUG -Wall -Wextra -pedantic 
#        -std=c++11 -pthread)
