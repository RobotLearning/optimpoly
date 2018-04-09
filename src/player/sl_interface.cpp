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
using namespace player;

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
 *
 */
struct blob_state {
	int status; //!< was ball detected reliably in cameras
	double pos[NCART]; //!< ball center cartesian positions from cameras 1 and 2(after calibration)
};

player_flags flags; //!< global structure for setting Player options

/*
 *
 * Fusing the two blobs
 * If both blobs are valid last blob is preferred
 * Only updates if the blobs are valid, i.e. not obvious outliers
 *
 */
static bool fuse_blobs(const blob_state blobs[NBLOBS], vec3 & obs);

/*
 * Checks for the validity of blob ball data using obvious table tennis checks.
 * Returns TRUE if valid.
 *
 * Only checks for obvious outliers. Does not use uncertainty estimates
 * to assess validity so do not rely on it as the sole source of outlier detection!
 *
 *
 * @param blob Blob structure. Contains a status boolean
 * variable and cartesian coordinates (indices zero-to-two).
 * @param verbose If verbose is TRUE, then detecting obvious outliers will
 * print to standard output.
 * @return valid If ball is valid (status is TRUE and not an obvious outlier)
 * return true.
 */
static bool check_blob_validity(const blob_state blobs[NBLOBS],
                                bool verbose);

/*
 *
 * Saves actual ball data if save flag is set to TRUE
 * and the ball observations and estimated ball state one another
 */
static void save_ball_data(const blob_state blobs[NBLOBS],
                           const Player *robot,
                           const KF & filter,
                           std::ofstream & stream);

/*
 *  Set algorithm to initialize Player with.
 *  alg_num selects between three algorithms: VHP/FOCUSED/DEFENSIVE.
 */
static void set_algorithm(const int alg_num);

#include "sl_interface.h"

void set_algorithm(const int alg_num) {

	switch (alg_num) {
		case 0:
			std::cout << "Setting to FOCUSED player..." << std::endl;
			flags.alg = FOCUS;
			break;
		case 1:
			std::cout << "Setting to DEFENSIVE player..." << std::endl;
			flags.alg = DP;
			break;
		case 2:
			std::cout << "Setting to VHP player..." << std::endl;
			flags.alg = VHP;
			break;
		default:
			flags.alg = FOCUS;
	}
}

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
			("penalty_loc", po::value<std::vector<double>>(&flags.penalty_loc)->multitoken(), "punishment locations for DP")
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

static bool check_blob_validity(const blob_state & blob, bool verbose) {

	bool valid = true;
	static double last_blob[NCART];
	static double zMax = 0.5;
	static double zMin = floor_level - table_height;
	static double xMax = table_width/2.0;
	static double yMax = 0.5;
	static double yMin = dist_to_table - table_length - 1.0;
	static double yCenter = dist_to_table - table_length/2.0;

	// if both blobs are valid
	// then check for both
	// else check only one

	if (blob.status == false) {
		valid = false;
	}
	else if (blob.pos[Z] > zMax) {
		if (verbose)
			printf("BLOB NOT VALID! Ball is above zMax = 0.5!\n");
		valid = false;
	}
	// on the table blob should not appear under the table
	else if (fabs(blob.pos[X]) < xMax && fabs(blob.pos[Y] - yCenter) < table_length/2.0
			&& blob.pos[Z] < zMin) {
		if (verbose)
			printf("BLOB NOT VALID! Ball appears under the table!\n");
		valid = false;
	}
	else if (last_blob[Y] < yCenter && blob.pos[Y] > dist_to_table) {
		if (verbose)
			printf("BLOB NOT VALID! Ball suddenly jumped in Y!\n");
		valid = false;
	}
	else if (blob.pos[Y] < yMin || blob.pos[Y] > yMax) {
		if (verbose)
			printf("BLOB NOT VALID! Ball is outside y limits [%f,%f]!\n", yMin, yMax);
		valid = false;
	}

	last_blob[X] = blob.pos[X];
	last_blob[Y] = blob.pos[Y];
	last_blob[Z] = blob.pos[Z];
	return valid;
}

static bool fuse_blobs(const blob_state blobs[NBLOBS], vec3 & obs) {

	bool status = false;

	// if ball is detected reliably
	// Here we hope to avoid outliers and prefer the blob3 over blob1
	if (check_blob_validity(blobs[0],flags.verbosity > 2) ||
			check_blob_validity(blobs[1],flags.verbosity > 2)) {
		status = true;
		if (blobs[1].status) { //cameras 3 and 4
			for (int i = X; i <= Z; i++)
				obs(i) = blobs[1].pos[i];
		}
		else { // cameras 1 and 2
			for (int i = X; i <= Z; i++)
				obs(i) = blobs[0].pos[i];
		}
	}
	return status;
}

void play(const SL_Jstate joint_state[NDOF+1],
		  const blob_state blobs[NBLOBS],
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

static void save_ball_data(const blob_state blobs[NBLOBS],
                            const Player *robot,
                            const KF & filter,
                            std::ofstream & stream) {

	static rowvec ball_full;
	vec6 ball_est = zeros<vec>(6);

	if (flags.save && (blobs[0].status || blobs[1].status)) {

		if (robot->filter_is_initialized()) {
			ball_est = filter.get_mean();
		}
		ball_full << 1 << blobs[0].status << blobs[0].pos[X] << blobs[0].pos[Y] << blobs[0].pos[Z]
				  << 3 << blobs[1].status << blobs[1].pos[X] << blobs[1].pos[Y] << blobs[1].pos[Z] << endr;
		//cout << ball_full << endl;
		ball_full = join_horiz(ball_full,ball_est.t());
		if (stream.is_open()) {
			stream << ball_full;
		}
		//stream_balls.close();
	}
}

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
