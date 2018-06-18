cmake_minimum_required (VERSION 2.8)
project(${PROJECT_NAME})

find_library(ZMQ zmq)
# Using BOOST UNIT TEST FRAMEWORK for testing

# CREATE SHARED LIBRARY
set(TEST_EXEC unit_tests)
set(TEST_SRC
    test_kalman.cpp
    test_kinematics.cpp
    test_optim.cpp
    test_table_tennis.cpp
    test_interface.cpp
    test_serve.cpp
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

find_library(UNIT_TEST NAMES libboost_unit_test_framework.a)
#MESSAGE(STATUS "BOOST UNIT TEST LOCATION: " ${UNIT_TEST})
# INCLUDE "OUR" SHARED LIBRARY
target_link_libraries(${TEST_EXEC}
    ${UNIT_TEST}
    ${ZMQPP}
    ${ZMQ}
    ${PROJECT_NAME}
    armadillo
    boost_program_options)

# INSTALL FOLDER
install(TARGETS ${TEST_EXEC}
    DESTINATION ${CMAKE_SOURCE_DIR})

# OPTIMIZATION OPTIONS
#target_compile_options(${TEST_EXEC} PRIVATE
#                -O3 -Wall -DNDEBUG -std=c++11 -pthread)

# DEBUG OPTIONS
#target_compile_options(example PUBLIC 
#        -g -DDEBUG -Wall -Wextra -pedantic 
#        -std=c++11 -pthread)
