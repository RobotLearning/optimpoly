cmake_minimum_required(VERSION 2.8)

set(TASKNAME table_tennis)

project(${TASKNAME})

message( STATUS "Adding ${TASKNAME} task for ${NAME}" )

# add the following line if the task requires to communicate with Matlab
#required_matlab_shared_mem()


include_directories(
  #${CMAKE_SOURCE_DIR}/IAS/include
  ~/table-tennis/include
  ~/table-tennis/include/player
  ~/table-tennis/include/optim
)

# for executables
set(ADD_SRCS_TASK
    #EXAMPLE_sim.c
    table_tennis_task.c
    table_tennis_common.c
    # any other additional file needed
)

# for graphics
set(ADD_SRCS_OPENGL
    #table_tennis_graphics.c
    # any other additional file needed
)

link_directories(~/table-tennis/lib)
link_libraries(~/table-tennis/lib/libtable-tennis.so nlopt)

finalize_auto_add_task()
