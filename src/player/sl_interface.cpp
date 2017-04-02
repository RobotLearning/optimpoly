/**
 *
 * @file sl_interface.cpp
 *
 * @brief Interface of the Player class to the SL real-time simulator and to
 * the robot.
 *
 */

#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <sys/time.h>
#include "kalman.h"
#include "player.hpp"
#include "tabletennis.h"

using namespace arma;

/* The data structures from SL */
/**
 * @brief (actual) joint space state for each DOF
 */
struct SL_Jstate {
	double   th;   /*!< theta */
	double   thd;  /*!< theta-dot */
	double   thdd; /*!< theta-dot-dot */
	double   ufb;  /*!< feedback portion of command */
	double   u;    /*!< torque command */
	double   load; /*!< sensed torque */
};
/**
 * @brief (desired) joint space state commands for each DOF
 */
struct SL_DJstate { /*!< desired values for controller */
	double   th;   /*!< theta */
	double   thd;  /*!< theta-dot */
	double   thdd; /*!< theta-dot-dot */
	double   uff;  /*!< feedforward torque command */
	double   uex;  /*!< externally imposed torque */
};

/**
 * @brief (actual) Cartesian state
 */
struct SL_Cstate {
	double   x[NCART+1];    /*!< Position [x,y,z] */
	double   xd[NCART+1];   /*!< Velocity */
	double   xdd[NCART+1];  /*!< Acceleration */
};

/**
 * @brief Vision blob coming from SL.
 */
struct SL_VisionBlob {
	char       status;
	SL_Cstate  blob;
};

/**
 * @brief Options passed to Player class (algorithm, saving, corrections).
 */
struct pflags { //! player flags
	int verbosity = 0; //! OFF, LOW, HIGH
	bool reset = true; //! reinitializing player class
	bool save = false; //! saving ball/robot data
	bool mpc = false; //! turn on/off corrections
	algo alg = FIXED; //! algorithm for trajectory generation
};

pflags player_flags; //! global structure for setting Player options

#include "sl_interface.h"

/**
 * @brief Set algorithm and options to initialize Player with.
 *
 * The global variable player_flags is set here and
 * the play() function will use it to initialize the Player class.
 *
 * @param alg_num Select between three algorithms: VHP/FIXED/LAZY.
 * @param mpc_flag Turn on/off online corrections.
 * @param save_flag Save ball observations if TRUE.
 * @param verbose_flag Flag for (verbose) printing, 0 = OFF, 1 = LOW, 2 = HIGH.
 */
void set_algorithm(const int alg_num, const int mpc_flag,
		           const int save_flag, const int verbose_flag) {

	player_flags.reset = true;
	switch (alg_num) {
		case 0:
			std::cout << "Setting to VHP player..." << std::endl;
			player_flags.alg = VHP;
			break;
		case 1:
			std::cout << "Setting to FIXED player..." << std::endl;
			player_flags.alg = FIXED;
			break;
		case 2:
			std::cout << "Setting to LAZY player..." << std::endl;
			player_flags.alg = LAZY;
			break;
		default:
			player_flags.alg = FIXED;
	}
	mpc_flag ? player_flags.mpc = true : player_flags.mpc = false;
	save_flag ? player_flags.save = true : player_flags.save = false;
	player_flags.verbosity = verbose_flag;

}

/*
 *
 * Checks for the validity of blob ball data using obvious table tennis checks.
 * Returns TRUE if valid.
 *
 * Only checks for obvious outliers. Does not use uncertainty estimates
 * to assess validity so do not rely on it as the sole source of outlier detection!
 *
 *
 * @param blob Blob structure from SL (one-indexed). Contains a status boolean
 * variable and cartesian coordinates (indices one-to-three).
 * @param verbose If verbose is TRUE, then detecting obvious outliers will
 * print to standard output.
 * @return valid If ball is valid (status is TRUE and not an obvious outlier)
 * return true.
 */
static bool check_blob_validity(const SL_VisionBlob & blob, bool verbose) {

	bool valid = true;
	static double zMax = 0.5;
	static double zMin = floor_level - table_height;
	static double xMax = table_width/2.0;
	static double yMax = 0.5;
	static double yMin = dist_to_table - table_length - 1.0;
	static double yCenter = dist_to_table - table_length/2.0;

	if (blob.status == false) {
		if (verbose)
			printf("BLOB NOT VALID! Ball status is false!\n");
		valid = false;
	}
	else if (blob.blob.x[3] > zMax) {
		if (verbose)
			printf("BLOB NOT VALID! Ball is above zMax = 0.5!\n");
		valid = false;
	}
	else if (blob.blob.x[2] > yMax || blob.blob.x[2] < yMin) {
		if (verbose)
			printf("BLOB NOT VALID! Ball is outside y-limits!\n");
		valid = false;
	}
	// between the table if ball appears outside the x limits
	else if (fabs(blob.blob.x[2] - yCenter) < table_length/2.0 &&
			fabs(blob.blob.x[1]) > xMax) {
		if (verbose)
			printf("BLOB NOT VALID! Ball is outside the table!\n");
		valid = false;
	}
	// on the table blob should not appear under the table
	else if (fabs(blob.blob.x[1]) < xMax && fabs(blob.blob.x[2] - yCenter) < table_length/2.0
			&& blob.blob.x[3] < zMin) {
		if (verbose)
			printf("BLOB NOT VALID! Ball appears under the table!\n");
		valid = false;
	}
	return valid;
}

/*
 *
 * Fusing the blobs
 * If both blobs are valid blob3 is preferred
 * Only updates if the blobs are valid, i.e. not obvious outliers
 *
 */
static bool fuse_blobs(const SL_VisionBlob blobs[4], vec3 & obs) {

	static bool status = false;

	// if ball is detected reliably
	// Here we hope to avoid outliers and prefer the blob3 over blob1
	if (check_blob_validity(blobs[3],player_flags.verbosity) ||
			check_blob_validity(blobs[1],player_flags.verbosity)) {
		status = true;
		if (blobs[3].status) {
			for (int i = X; i <= Z; i++)
				obs(i) = blobs[3].blob.x[i+1];
		}
		else {
			for (int i = X; i <= Z; i++)
				obs(i) = blobs[1].blob.x[i+1];
		}
	}
	return status;
}

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
void play(const SL_Jstate joint_state[NDOF+1],
		  const SL_VisionBlob blobs[4],
		  SL_DJstate joint_des_state[NDOF+1]) {

	static vec7 q0;
	static vec3 ball_obs;
	static joint qact;
	static joint qdes;
	static Player *robot = nullptr; // centered player
	static EKF filter = init_filter();

	if (player_flags.reset) {
		for (int i = 0; i < NDOF; i++) {
			qdes.q(i) = q0(i) = joint_state[i+1].th;
			qdes.qd(i) = 0.0;
			qdes.qdd(i) = 0.0;
		}
		filter = init_filter();
		delete robot;
		robot = new Player(q0,filter,player_flags.alg,player_flags.mpc,player_flags.verbosity);
		player_flags.reset = false;
	}
	else {
		for (int i = 0; i < NDOF; i++) {
			qact.q(i) = joint_state[i+1].th;
			qact.qd(i) = joint_state[i+1].thd;
			qact.qdd(i) = joint_state[i+1].thdd;
		}
		fuse_blobs(blobs,ball_obs);
		robot->play(qact,ball_obs,qdes);
		save_data(qact,qdes,blobs,ball_obs,filter);
	}

	// update desired joint state
	for (int i = 0; i < NDOF; i++) {
		joint_des_state[i+1].th = qdes.q(i);
		joint_des_state[i+1].thd = qdes.qd(i);
		joint_des_state[i+1].thdd = qdes.qdd(i);
	}

}

/*
 *
 * Saves actual and desired joints to one file if save flag is set to TRUE
 * and the ball observations and estimated ball state one another
 *
 * TODO: no need to open close each time!
 *
 */
static void save_data(const joint & qact, const joint & qdes,
		       const SL_VisionBlob blobs[4], const vec3 & ball_obs, const KF & filter) {

	static rowvec ball_full;
	static std::string home = std::getenv("HOME");
	static std::string joint_file = home + "/polyoptim/joints.txt";
	static std::string ball_file = home + "/polyoptim/balls.txt";
	static std::ofstream stream_joints;
	static std::ofstream stream_balls;
	static vec6 ball_est = zeros<vec>(6);

	if (player_flags.save) {
		try {
			ball_est = filter.get_mean();
		}
		catch (const char * exception) {
			// do nothing
		}

		/*rowvec qdes_full = join_horiz(join_horiz(qdes.q.t(),qdes.qd.t()),qdes.qdd.t());
		rowvec qact_full = join_horiz(join_horiz(qact.q.t(),qact.qd.t()),qact.qdd.t());
		rowvec q_full = join_horiz(qdes_full, qact_full);
		stream_joints.open(joint_file,ios::out | ios::app);
		if (stream_joints.is_open()) {
			stream_joints << q_full;
		}*/

		stream_balls.open(ball_file,ios::out | ios::app);
		ball_full << 1 << ((int)blobs[1].status)
				  << blobs[1].blob.x[1] << blobs[1].blob.x[2] << blobs[1].blob.x[3]
				  << 3 << ((int)blobs[3].status)
				  << blobs[3].blob.x[1] << blobs[3].blob.x[2] << blobs[3].blob.x[3] << endr;
		ball_full = join_horiz(join_horiz(ball_full,ball_obs.t()),ball_est.t());
		if (stream_balls.is_open()) {
			stream_balls << ball_full;
		}
		//stream_joints.close();
		stream_balls.close();
	}
}

/**
 * @brief  CHEAT with exact knowledge of ball state.
 *
 * Interface to the PLAYER class that generates desired hitting trajectories.
 * First initializes the player and then starts calling cheat() interface function.
 *
 * @param joint_state Actual joint positions, velocities, accelerations.
 * @param sim_ball_state Exact simulated ball state (positions and velocities).
 * @param joint_des_state Desired joint position, velocity and acceleration commands.
 */
void cheat(const SL_Jstate joint_state[NDOF+1],
		  const SL_Cstate sim_ball_state,
		  SL_DJstate joint_des_state[NDOF+1]) {

	static vec7 q0;
	static vec6 ball_state;
	static joint qact;
	static joint qdes;
	static Player *cp; // centered player
	static EKF filter = init_filter();

	if (player_flags.reset) {
		for (int i = 0; i < NDOF; i++) {
			qdes.q(i) = q0(i) = joint_state[i+1].th;
			qdes.qd(i) = 0.0;
			qdes.qdd(i) = 0.0;
		}
		cp = new Player(q0,filter,player_flags.alg,player_flags.mpc,player_flags.verbosity);
		player_flags.reset = false;
	}
	else {
		for (int i = 0; i < NDOF; i++) {
			qact.q(i) = joint_state[i+1].th;
			qact.qd(i) = joint_state[i+1].thd;
			qact.qdd(i) = joint_state[i+1].thdd;
		}
		for (int i = 0; i < NCART; i++) {
			ball_state(i) = sim_ball_state.x[i+1];
			ball_state(i+NCART) = sim_ball_state.xd[i+1];
		}
		cp->cheat(qact,ball_state,qdes);
	}

	// update desired joint state
	for (int i = 0; i < NDOF; i++) {
		joint_des_state[i+1].th = qdes.q(i);
		joint_des_state[i+1].thd = qdes.qd(i);
		joint_des_state[i+1].thdd = qdes.qdd(i);
	}
}


// * Inverts the given SL matrix matc [with indices starting from 1]
// * by initializing ARMADILLO equivalent and inverting it
// *
// * After inversion the matrix contents are overwritten to the SL input matrix (double**)
// * Memory for output matrix needs to be allocated before in SL (e.g. my_matrix()).
// *
// */
//void invert_matrix(double** matc, int nrows, double** out) {
//
//	double* bowels = new double[nrows*nrows];
//
//	// fill the bowels
//	int i,j;
//	for(i = 1; i <= nrows; i++) {
//		for(j = 1; j <= nrows; j++) {
//			bowels[(j-1)*nrows + (i-1)] = matc[i][j]; // mat from SL
//		}
//	}
//
//	// initialize ARMA mat
//	mat A(bowels,nrows,nrows);
//	mat B = inv(A);
//
//	for(i = 1; i <= nrows; i++) {
//		for(j = 1; j <= nrows; j++) {
//			out[i][j] = B(i-1,j-1); // mat from ARMADILLO
//		}
//	}
//
//}
//

// * Taking pseudoinverse with default tolerance using ARMADILLO
// *
// * Output matrix must be the transpose of the input matrix (so make sure
// * you initialize properly!).
// *
// *
// */
//void pinv_matrix(double** matc, int nrows, int ncols, double** out) {
//
//	double* bowels = new double[nrows*ncols];
//
//	// fill the bowels
//	int i,j;
//	for(i = 1; i <= nrows; i++) {
//		for(j = 1; j <= ncols; j++) {
//			bowels[(j-1)*nrows + (i-1)] = matc[i][j]; // mat from SL
//		}
//	}
//
//	// initialize ARMA mat
//	mat A(bowels,nrows,ncols);
//	mat B = pinv(A);
//
//	for(i = 1; i <= ncols; i++) {
//		for(j = 1; j <= nrows; j++) {
//			out[i][j] = B(i-1,j-1); // mat from ARMADILLO
//		}
//	}
//}

