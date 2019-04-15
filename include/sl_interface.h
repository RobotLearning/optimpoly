/**
 * @file sl_interface.h
 * @brief Interface to SL for both play and serve tasks.
 * Player Exposes play() and cheat() modes.
 *
 * Before starting play() mode, relevant options for the Player class
 * must be set by calling load_player_options().
 *
 */
#pragma once
#ifdef __cplusplus

#include "sl_structs.h"
// nothing here

#endif

#ifdef __cplusplus
extern "C" {
#endif

/** @brief New interface to PLAYER class. Uses ZMQ connection to fetch 3d ball info.
 *
 * Generates desired hitting trajectories.
 *
 * First initializes the player according to the pre-set options
 * and then starts calling play() interface function. Must be called every DT ms.
 *
 * @param joint_state Actual joint positions, velocities, accelerations.
 * @param joint_des_state Desired joint position, velocity and acceleration commands. */
extern void play(const SL_Jstate joint_state[],
                 SL_DJstate joint_des_state[]);


/**
 * @brief CHEAT with exact knowledge of ball state.
 *
 * Interface to the PLAYER class that generates desired hitting trajectories.
 * First initializes the player and then starts calling cheat() interface function.
 *
 * @param joint_state Actual joint positions, velocities, accelerations.
 * @param sim_ball_state Exact simulated ball state (positions and velocities).
 * @param joint_des_state Desired joint position, velocity and acceleration commands.
 */
extern void cheat(const SL_Jstate joint_state[],
                  const SL_Cstate sim_ball_state,
                  SL_DJstate joint_des_state[]);

/**
 * @brief Set algorithm and options to initialize Player with.
 *
 * The global variable flags is set here and
 * the play() function will use it to initialize the Player class.
 *
 */
extern void* load_player_options();

/**
 * @brief Called from SL, save joint data to a file. If save_qdes is true, then save also the desired joints.
 * If reset is true, then close and open the stream.
 * */
extern void save_joint_data(const SL_Jstate joint_state[],
                     const SL_DJstate joint_des_state[],
                     const int save_qdes,
					 const int reset);

/** @brief Called from SL, creates a Listener with given URL and saves data that it can fetch with it */
extern void save_ball_data(const char* url_string,
							const int debug_vision,
							const int listen_2d,
							const int reset);

/* SERVE FUNCTIONS AND STRUCTS EXPOSED */

typedef struct {
    int use_inv_dyn_fb; //!< use computed-torque if false
} serve_task_options;

/** \brief Set options for the serve, i.e. ServeBall class */
extern void* load_serve_options(double custom_pose[], serve_task_options *opt);

/** \brief Serve with learned DMPs */
extern void serve_ball(const SL_Jstate joint_state[],
                     SL_DJstate joint_des_state[]);


#ifdef __cplusplus
} // extern "C"
#endif
