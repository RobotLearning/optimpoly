cmake_minimum_required (VERSION 2.8)
project(${PROJECT_NAME})

# CREATE SHARED LIBRARY
set(SOURCES
    player/extkalman.cpp
    player/kalman.cpp
    player/kinematics.cpp
    player/lookup.cpp
    player/player.cpp
    player/sl_interface.cpp
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
)
add_library (${PROJECT_NAME} SHARED ${SOURCES})

# SET PROJECT VERSION
set_target_properties(${PROJECT_NAME} PROPERTIES 
    VERSION ${PROJECT_VERSION})
#    SOVERSION 1
#    PUBLIC_HEADER api/mylib.h)

# INCLUDE HEADERS (top folder, as opposed to CMAKE_CURRENT_SOURCE_DIR)
target_include_directories (${PROJECT_NAME} PRIVATE
    ${CMAKE_SOURCE_DIR}/include/optim
    ${CMAKE_SOURCE_DIR}/include/player)

# LINK EXTERNAL LIBRARIES HERE
target_link_libraries(${PROJECT_NAME} armadillo boost_program_options)

# INSTALL FOLDER
install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib)
#    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# OPTIMIZATION OPTIONS
#target_compile_options(${PROJECT_NAME} PRIVATE
#                -O3 -Wall -DNDEBUG -std=c++11 -pthread)

# DEBUG OPTIONS
#target_compile_options(example PUBLIC 
#        -g -DDEBUG -Wall -Wextra -pedantic 
#        -std=c++11 -pthread)
