cmake_minimum_required (VERSION 2.8)
project(${PROJECT_NAME})
find_library(ZMQ zmq)
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

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

target_link_libraries(${TEST_EXEC}
    ${GTEST_LIBRARY}
    ${ZMQPP}
    ${ZMQ}
    ${PROJECT_NAME}
    ${ARMADILLO_LIBRARY}
    pthread
)

# INSTALL FOLDER
install(TARGETS ${TEST_EXEC}
    DESTINATION ${CMAKE_SOURCE_DIR})
