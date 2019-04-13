cmake_minimum_required (VERSION 2.8)
project(${PROJECT_NAME})

# CREATE SHARED LIBRARY
set(SOURCES
    player/extkalman.cpp
    player/kalman.cpp
    player/kinematics.cpp
    player/lookup.cpp
    player/player.cpp
    player/table_tennis.cpp
    player/traj.cpp
    optim/defensive_optim.cpp
    optim/estimate_ball.cpp
    optim/focused_optim.cpp
    optim/kinematics.cpp
    optim/optim.cpp
    optim/racket_optim.cpp
    optim/rest_optim.cpp
    optim/utils.cpp
    optim/vhp_optim.cpp
    serve/dmp.cpp
    serve/serve.cpp
    serve/rbf.cpp
    sl_interface.cpp
    ball_interface.cpp
)
add_library (${PROJECT_NAME} SHARED ${SOURCES})

# SET PROJECT VERSION
set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION ${PROJECT_VERSION})

# INCLUDE HEADERS (top folder, as opposed to CMAKE_CURRENT_SOURCE_DIR)
target_include_directories (${PROJECT_NAME} PRIVATE
    ${ZMQPP_INCLUDES}
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/serve
    ${CMAKE_SOURCE_DIR}/include/optim
    ${CMAKE_SOURCE_DIR}/include/player)

#message(STATUS "armadillo lib path = " ${ARMADILLO_LIBRARY})
#message(STATUS "boost program options lib path = " ${Boost_PROGRAM_OPTIONS_LIBRARY})

# LINK EXTERNAL LIBRARIES HERE
target_link_libraries(${PROJECT_NAME}
  ${ARMADILLO_LIBRARY}
  nlopt
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${ZMQPP})

# INSTALL FOLDER
install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib)
