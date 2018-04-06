/**
 * @file sl_interface.h
 * @brief Interface to SL. Exposes play() and cheat() modes.
 *
 * Before starting play() mode, relevant options for the Player class
 * must be set by calling set_algorithm().
 *
 */
#ifndef SL_INTERF_H
#define SL_INTERF_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Interface to the PLAYER class that generates desired hitting trajectories.
 *
 * First initializes the player according to the pre-set options
 * and then starts calling play() interface function. Must be called every DT ms.
 *
 *
 * @param joint_state Actual joint positions, velocities, accelerations.
 * @param blobs Two ball 3d-positions from 4-cameras are stored in blobs[1] and blobs[3]
 * @param joint_des_state Desired joint position, velocity and acceleration commands.
 */
extern void play(const SL_Jstate joint_state[],
		         const blob_state blobs[],
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
extern void load_options();

#ifdef __cplusplus
} // extern "C"
#endif

#ifdef __cplusplus

#endif

#endif /* SL_INTERF_H */
