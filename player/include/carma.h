
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

extern void set_algorithm(int num);

// inverts a given matrix and overwriting contents with inverse
extern void invert_matrix(double** mat, int nrows, double** out);
// pseudoinverse of given matrix matc
extern void pinv_matrix(double** matc, int nrows, int ncols, double** out);


#ifdef __cplusplus
} // extern "C"
#endif

#ifdef __cplusplus
// internal c++ functions
static bool fuse_blobs(const SL_VisionBlob blobs[], vec3 & obs);
static bool check_blob_validity(const SL_VisionBlob & blob, bool verbose);

#endif

#endif /* CARMA_H */
