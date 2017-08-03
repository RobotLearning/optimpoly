/**
 *
 * @file sl_interface.cpp
 *
 * @brief Interface of the Player class to the SL real-time simulator and to
 * the robot.
 *
 */

#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
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
 * @brief Vision blob info coming from SL (after calibration).
 */
struct SL_VisionBlob {
	char       status; //!< was ball detected reliably?
	SL_Cstate  blob; //!< ball center cartesian positions (after calibration)
};

player_flags flags; //!< global structure for setting Player options

#include "sl_interface.h"

/**
 * @brief Set algorithm to initialize Player with.
 *
 * @param alg_num Select between three algorithms: VHP/FIXED/LAZY.
 */
void set_algorithm(const int alg_num) {

	switch (alg_num) {
		case 0:
			std::cout << "Setting to FOCUSED player..." << std::endl;
			flags.alg = FOCUS;
			break;
		case 1:
			std::cout << "Setting to LAZY player..." << std::endl;
			flags.alg = LAZY;
			break;
		case 2:
			std::cout << "Setting to VHP player..." << std::endl;
			flags.alg = VHP;
			break;
		default:
			flags.alg = FOCUS;
	}
}

// A helper function to simplify the main part.
//template<class T>
//std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
//    copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
//    return os;
//}

/**
 * @brief Set algorithm and options to initialize Player with.
 *
 * The global variable flags is set here and
 * the play() function will use it to initialize the Player class.
 *
 */
void load_options() {

	namespace po = boost::program_options;
	using namespace std;

	flags.reset = true;
	string home = std::getenv("HOME");
	string config_file = home + "/polyoptim/" + "player.cfg";
	int alg_num;

    try {
		// Declare a group of options that will be
		// allowed in config file
		po::options_description config("Configuration");
		config.add_options()
		    ("outlier_detection", po::value<bool>(&flags.outlier_detection)->default_value(true),
			      "OUTLIER DETECTION FOR REAL ROBOT!")
		    ("rejection_multiplier", po::value<double>(&flags.out_reject_mult)->default_value(2.0),
				"OUTLIER DETECTION MULTIPLIER FOR REAL ROBOT!")
			("check_bounce", po::value<bool>(&flags.check_bounce),
			      "checking bounce before moving!")
		    ("weights", po::value<std::vector<double>>(&flags.weights)->multitoken(), "hit,net,land weights for DP")
			("mult_vel", po::value<std::vector<double>>(&flags.mult_vel)->multitoken(), "velocity mult. for DP")
			("rest_posture_optim", po::value<bool>(&flags.optim_rest_posture)->default_value(false),
				"turn on resting state optimization")
			//("lookup", po::value<bool>(&flags.lookup), "start moving with lookup")
			("algorithm", po::value<int>(&alg_num)->default_value(0),
				  "optimization method")
			("mpc", po::value<bool>(&flags.mpc)->default_value(false),
				 "corrections (MPC)")
			("spin", po::value<bool>(&flags.spin)->default_value(false),
						 "apply spin model")
			("verbose", po::value<int>(&flags.verbosity)->default_value(1),
		         "verbosity level")
		    ("save_data", po::value<bool>(&flags.save)->default_value(false),
		         "saving robot/ball data")
		    ("ball_land_des_x_offset", po::value<double>(&flags.ball_land_des_offset[0]),
		    	 "ball land x offset")
			("ball_land_des_y_offset", po::value<double>(&flags.ball_land_des_offset[1]),
				 "ball land y offset")
		    ("time_land_des", po::value<double>(&flags.time_land_des),
		    	 "time land des")
			("start_optim_offset", po::value<double>(&flags.optim_offset),
				 "start optim offset")
			("time2return", po::value<double>(&flags.time2return),
						 "time to return to start posture")
			("freq_mpc", po::value<int>(&flags.freq_mpc), "frequency of updates")
		    ("min_obs", po::value<int>(&flags.min_obs), "minimum obs to start filter")
		    ("var_noise", po::value<double>(&flags.var_noise), "std of filter obs noise")
		    ("var_model", po::value<double>(&flags.var_model), "std of filter process noise")
		    ("t_reset_threshold", po::value<double>(&flags.t_reset_thresh), "filter reset threshold time")
		    ("VHPY", po::value<double>(&flags.VHPY), "location of VHP");
        po::variables_map vm;
        ifstream ifs(config_file.c_str());
        if (!ifs) {
            cout << "can not open config file: " << config_file << "\n";
        }
        else {
            po::store(parse_config_file(ifs, config), vm);
            notify(vm);
        }
    }
    catch(exception& e) {
        cout << e.what() << "\n";
    }
    set_algorithm(alg_num);
    flags.detach = true; // always detached in SL/REAL ROBOT!
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
	static double last_blob[NCART+1];
	static double zMax = 0.5;
	static double zMin = floor_level - table_height;
	static double xMax = table_width/2.0;
	static double yMax = 0.5;
	static double yMin = dist_to_table - table_length - 1.0;
	static double yCenter = dist_to_table - table_length/2.0;

	if (blob.status == false) {
		//if (verbose)
		//	printf("BLOB NOT VALID! Ball status is false!\n");
		valid = false;
	}
	else if (blob.blob.x[3] > zMax) {
		if (verbose)
			printf("BLOB NOT VALID! Ball is above zMax = 0.5!\n");
		valid = false;
	}
	else if (last_blob[2] < yCenter && blob.blob.x[2] > dist_to_table) {
		if (verbose)
			printf("BLOB NOT VALID! Ball suddenly jumped in Y!\n");
		valid = false;
	}
	else if (blob.blob.x[2] < yMin || blob.blob.x[2] > yMax) {
		if (verbose)
			printf("BLOB NOT VALID! Ball is outside y limits [%f,%f]!\n", yMin, yMax);
		valid = false;
	}
	// on the table blob should not appear under the table
	else if (fabs(blob.blob.x[1]) < xMax && fabs(blob.blob.x[2] - yCenter) < table_length/2.0
			&& blob.blob.x[3] < zMin) {
		if (verbose)
			printf("BLOB NOT VALID! Ball appears under the table!\n");
		valid = false;
	}
	last_blob[1] = blob.blob.x[1];
	last_blob[2] = blob.blob.x[2];
	last_blob[3] = blob.blob.x[3];
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
	if (check_blob_validity(blobs[3],flags.verbosity > 2) ||
			check_blob_validity(blobs[1],flags.verbosity > 2)) {
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
	static std::ofstream stream_balls;
	static std::string home = std::getenv("HOME");
	static std::string ball_file = home + "/polyoptim/balls.txt";
	static EKF filter = init_filter(0.3,0.001,flags.spin);
	static int firsttime = true;

	if (firsttime && flags.save) {
		stream_balls.open(ball_file,std::ofstream::out | std::ofstream::app);
		firsttime = false;
	}

	if (flags.reset) {
		for (int i = 0; i < NDOF; i++) {
			qdes.q(i) = q0(i) = joint_state[i+1].th;
			qdes.qd(i) = 0.0;
			qdes.qdd(i) = 0.0;
		}
		filter = init_filter(0.3,0.001,flags.spin);
		delete robot;
		robot = new Player(q0,filter,flags);
		flags.reset = false;
	}
	else {
		for (int i = 0; i < NDOF; i++) {
			qact.q(i) = joint_state[i+1].th;
			qact.qd(i) = joint_state[i+1].thd;
			qact.qdd(i) = joint_state[i+1].thdd;
		}
		fuse_blobs(blobs,ball_obs);
		robot->play(qact,ball_obs,qdes);
		save_ball_data(blobs,robot,filter,stream_balls);
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
 * Saves actual ball data if save flag is set to TRUE
 * and the ball observations and estimated ball state one another
 *
 *
 */
static void save_ball_data(const SL_VisionBlob blobs[4], const Player *robot, const KF & filter, std::ofstream & stream) {

	static rowvec ball_full;
	vec6 ball_est = zeros<vec>(6);
	int status1 = (int)blobs[1].status;
	int status3 = (int)blobs[3].status;

	if (flags.save && (status1 || status3)) {

		if (robot->filter_is_initialized()) {
			ball_est = filter.get_mean();
		}
		ball_full << 1 << status1  << blobs[1].blob.x[1] << blobs[1].blob.x[2] << blobs[1].blob.x[3]
				  << 3 << status3  << blobs[3].blob.x[1] << blobs[3].blob.x[2] << blobs[3].blob.x[3] << endr;
		//cout << ball_full << endl;
		ball_full = join_horiz(ball_full,ball_est.t());
		if (stream.is_open()) {
			stream << ball_full;
		}
		//stream_balls.close();
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

	if (flags.reset) {
		for (int i = 0; i < NDOF; i++) {
			qdes.q(i) = q0(i) = joint_state[i+1].th;
			qdes.qd(i) = 0.0;
			qdes.qdd(i) = 0.0;
		}
		cp = new Player(q0,filter,flags);
		flags.reset = false;
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
