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

// Interface for Player
extern void play(const SL_Jstate joint_state[],
		         const blob_state blobs[],
				 SL_DJstate joint_des_state[]);

extern void cheat(const SL_Jstate joint_state[],
		  const SL_Cstate sim_ball_state,
		  SL_DJstate joint_des_state[]);

extern void load_options();

#ifdef __cplusplus
} // extern "C"
#endif

#ifdef __cplusplus
// internal c++ functions
static bool fuse_blobs(const blob_state blobs[NBLOBS], vec3 & obs);
static bool check_blob_validity(const blob_state blobs[NBLOBS], bool verbose);
static void save_ball_data(const blob_state blobs[NBLOBS], const Player *robot, const KF & filter, std::ofstream & stream);
static void set_algorithm(const int alg_num);

#endif

#endif /* SL_INTERF_H */
