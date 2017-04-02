/**
 * @file sl_interface.h
 * @brief Interface to SL. Exposes play() and cheat() modes.
 *
 * Before starting play() mode, relevant options for the Player class
 * must be set by calling set_algorithm().
 *
 */
#ifndef CARMA_H
#define CARMA_H

#ifdef __cplusplus
extern "C" {
#endif

// Interface for Player
extern void play(const SL_Jstate joint_state[],
		         const SL_VisionBlob blobs[],
				 SL_DJstate joint_des_state[]);

extern void cheat(const SL_Jstate joint_state[],
		  const SL_Cstate sim_ball_state,
		  SL_DJstate joint_des_state[]);

extern void set_algorithm(const int alg_num, const int mpc_flag,
		                  const int save_flag, const int verbosity);

#ifdef __cplusplus
} // extern "C"
#endif

#ifdef __cplusplus
// internal c++ functions
static bool fuse_blobs(const SL_VisionBlob blobs[], vec3 & obs);
static bool check_blob_validity(const SL_VisionBlob & blob, bool verbose);
static void save_data(const joint & qact, const joint & qdes,
		       const SL_VisionBlob blobs[4], const vec3 & ball_obs, const KF & filter);

#endif

#endif /* CARMA_H */
