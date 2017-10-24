cmake_minimum_required(VERSION 2.8)

set(TASKNAME player)

project(${TASKNAME})

message( STATUS "Adding ${TASKNAME} task for ${NAME}" )

# add the following line if the task requires to communicate with Matlab
#required_matlab_shared_mem()


include_directories(
	#${CMAKE_SOURCE_DIR}/IAS/include
	~/polyoptim/include/player
	~/polyoptim/include/optim
)

# for executables
set(ADD_SRCS_TASK
    #EXAMPLE_sim.c
    player_task.c
    table_tennis_common.c
    # any other additional file needed
)

# for graphics
set(ADD_SRCS_OPENGL
    #table_tennis_graphics.c
    # any other additional file needed
)

link_directories(~/polyoptim/lib)
link_libraries(~/polyoptim/lib/libplayer.so nlopt)

finalize_auto_add_task()
