/**
 * @file sl_interface.h
 * @brief Interface to SL. Exposes play() and cheat() modes.
 *
 * Before starting play() mode, relevant options for the Player class
 * must be set by calling load_options().
 *
 */
#ifndef SL_INTERF_H
#define SL_INTERF_H

#ifdef __cplusplus

#include "constants.h"
#include "sl_structs.h"
#include "ball_interface.h"

struct serve_flags {
    bool save_joint_act_data = false;
    bool save_joint_des_data = false;
    bool use_inv_dyn_fb = false;
    double Tmax = 1.0; //!< time to evolve dmp if tau = 1.0
    std::string json_file = "dmp4.json"; //!< json file to load dmp from
};

#endif

#ifdef __cplusplus
extern "C" {
#endif

/** @brief New interface to PLAYER class. Uses ZMQ connection to fetch 3d ball info.*/
extern void play_new(const SL_Jstate joint_state[],
                 SL_DJstate joint_des_state[]);

/**
 * @brief Old interface to the PLAYER class that generates desired hitting trajectories.
 *
 * First initializes the player according to the pre-set options
 * and then starts calling play() interface function. Must be called every DT ms.
 *
 * Old interface written for the OLD setup. Time stamp is not provided.
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

/** \brief Save joint data to a file. If save_qdes is true, then save also the desired joints. */
extern void save_joint_data(const SL_Jstate joint_state[],
                     const SL_DJstate joint_des_state[],
                     const int save_qdes);

typedef struct {
    int init_dmp; //!< initialize dmp code if true
    int use_inv_dyn_fb; //!< use computed-torque if false
} dmp_task_options;

/** \brief Initialize DMP for the serve */
extern void init_dmp_serve(double custom_pose[], dmp_task_options *opt);

/** \brief Serve with learned DMPs */
extern void serve_with_dmp(const SL_Jstate joint_state[],
                            SL_DJstate joint_des_state[],
                            dmp_task_options *opt);


#ifdef __cplusplus
} // extern "C"
#endif

#endif /* SL_INTERF_H */
