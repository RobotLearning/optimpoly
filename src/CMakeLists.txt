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
    sl_interface.cpp
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

# LINK EXTERNAL LIBRARIES HERE
#add_library(boost_program_options SHARED IMPORTED)
#set_target_properties(boost_program_options PROPERTIES IMPORTED_LOCATION /home/robolab/software/lib/libboost_program_options.so)
target_link_libraries(${PROJECT_NAME}
  armadillo
  nlopt
  boost_program_options
  ${ZMQPP})

# INSTALL FOLDER
install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib)
