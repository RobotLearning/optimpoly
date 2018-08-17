#define NCART 3
#include "sl_structs.h"
#include "optim.h"
#include "player.hpp"
#include "kalman.h"
#include "constants.h"
#include "table.h"
#include "stdio.h"
#include "math.h"
#include <armadillo>

using namespace const_tt;
using namespace arma;
using namespace player;

player_flags pflags; //!< global structure for setting Player task options

/** @brief Vision ball info coming from SL (after calibration). */
struct blob_state {
    int status = 0; //!< was ball detected reliably in cameras and is it new data?
    double pos[3] = {0.0, 0.0, 0.0}; //!< ball center cartesian positions from two cameras
};

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
static bool check_blob_validity(const blob_state & blob, bool verbose);

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
    bool verbose = false;

    // if ball is detected reliably
    // Here we hope to avoid outliers and prefer the blob3 over blob1
    if (check_blob_validity(blobs[1],verbose) ||
            check_blob_validity(blobs[0],verbose)) {
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

void save_ball_data(const ball_obs & blobs,
                    const KF & filter) {

    static rowvec ball_full;
    static bool firsttime = true;
    static std::ofstream stream_balls;
    static std::string home = std::getenv("HOME");
    static std::string ball_file = home + "/table-tennis/balls.txt";
    if (firsttime) {
        stream_balls.open(ball_file,std::ofstream::out);
        firsttime = false;
    }
    vec6 ball_est = zeros<vec>(6);

    if (blobs.status) {

        try {
            ball_est = filter.get_mean();
            //cout << ball_full << endl;
            ball_full = join_horiz(blobs.pos,ball_est.t());
            if (stream_balls.is_open()) {
                stream_balls << ball_full;
            }
        }
        catch (const std::exception & not_init_error) {
            // filter not yet initialized
        }
        //stream_balls.close();
    }
}


void play(const SL_Jstate joint_state[NDOF+1],
		  const blob_state blobs[NBLOBS],
		  SL_DJstate joint_des_state[NDOF+1]) {

	static vec7 q0;
	static vec3 ball_obs;
	static optim::joint qact;
	static optim::joint qdes;
	static Player *robot = nullptr;
	static std::string home = std::getenv("HOME");
	static std::string ball_file = home + "/table-tennis/balls.txt";
	static EKF filter = init_ball_filter(0.3,0.001,pflags.spin);

	if (pflags.reset) {
		for (int i = 0; i < NDOF; i++) {
			qdes.q(i) = q0(i) = joint_state[i+1].th;
			qdes.qd(i) = 0.0;
			qdes.qdd(i) = 0.0;
		}
		filter = init_ball_filter(0.3,0.001,pflags.spin);
		delete robot;
		robot = new Player(q0,filter,pflags);
		pflags.reset = false;
	}
	else {
		for (int i = 0; i < NDOF; i++) {
			qact.q(i) = joint_state[i+1].th;
			qact.qd(i) = joint_state[i+1].thd;
			qact.qdd(i) = joint_state[i+1].thdd;
		}
		fuse_blobs(blobs,ball_obs);
		robot->play(qact,ball_obs,qdes);
		if (pflags.save)
		    save_ball_data(blobs,filter);
	}

	// update desired joint state
	for (int i = 0; i < NDOF; i++) {
		joint_des_state[i+1].th = qdes.q(i);
		joint_des_state[i+1].thd = qdes.qd(i);
		joint_des_state[i+1].thdd = qdes.qdd(i);
	}
}

void cheat(const SL_Jstate joint_state[NDOF+1],
          const SL_Cstate sim_ball_state,
          SL_DJstate joint_des_state[NDOF+1]) {

    static vec7 q0;
    static vec6 ball_state;
    static optim::joint qact;
    static optim::joint qdes;
    static Player *cp; // centered player
    static EKF filter = init_ball_filter();

    if (pflags.reset) {
        for (int i = 0; i < NDOF; i++) {
            qdes.q(i) = q0(i) = joint_state[i+1].th;
            qdes.qd(i) = 0.0;
            qdes.qdd(i) = 0.0;
        }
        cp = new Player(q0,filter,pflags);
        pflags.reset = false;
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
