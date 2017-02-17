
#ifndef CARMA_H
#define CARMA_H

#ifdef __cplusplus
extern "C" {
#endif

// kf function for SL
extern void ekf(double x[], double y[], double racket_pos[], double racket_orient[], int *reset);

// inverts a given matrix and overwriting contents with inverse
extern void invert_matrix(double** mat, int nrows, double** out);
// pseudoinverse of given matrix matc
extern void pinv_matrix(double** matc, int nrows, int ncols, double** out);


#ifdef __cplusplus
} // extern "C"
#endif

#ifdef __cplusplus
// internal c++ functions
void ekf(double x[], double y[], double racket_pos[], double racket_orient[], int *reset);
bool check_new_obs(const vec3 & obs);
EKF* init_ball_filter(const mat & obss, const vec & times);
void pass_mean_estimate(double x[], EKF * filter);
void estimate_prior(const mat & observations, const vec & times,
		            vec & x, mat & P);
#endif

#endif /* CARMA_H */
