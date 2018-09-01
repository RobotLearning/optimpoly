/**
 *
 * @file sl_interface.cpp
 *
 * @brief Interface of the Player class to the SL real-time simulator and to
 * the robot.
 *
 */

#include <boost/program_options.hpp>
#include "json.hpp"
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <armadillo>
#include <cmath>
#include <chrono>
#include <sys/time.h>
#include "kalman.h"
#include "player.hpp"
#include "dmp.h"
#include "rbf.h"
#include "serve.h"
#include "tabletennis.h"
#include "sl_interface.h"
#include "ball_interface.h"

using namespace arma;
using namespace player;
using namespace serve;
using namespace const_tt;

player_flags pflags; //!< global structure for setting Player task options
serve_flags sflags; //!< global structure for setting Serve task options

/*
 *  Set algorithm to initialize Player with.
 *  alg_num selects between three algorithms: VHP/FOCUSED/DEFENSIVE.
 */
static void set_player_algorithm(const int alg_num);

void set_player_algorithm(const int alg_num) {

	switch (alg_num) {
		case 0:
			std::cout << "Setting to FOCUSED player..." << std::endl;
			pflags.alg = FOCUS;
			break;
		case 1:
			std::cout << "Setting to DEFENSIVE player..." << std::endl;
			pflags.alg = DP;
			break;
		case 2:
			std::cout << "Setting to VHP player..." << std::endl;
			pflags.alg = VHP;
			break;
		default:
			pflags.alg = FOCUS;
	}
}

void* load_player_options() {

	namespace po = boost::program_options;
	using namespace std;

	pflags.reset = true;
	string home = std::getenv("HOME");
	string config_file = home + "/table-tennis/config/" + "player.cfg";
	int alg_num;

    try {
		// Declare a group of options that will be
		// allowed in config file
		po::options_description config("Configuration");
		config.add_options()
		    ("outlier_detection", po::value<bool>(&pflags.outlier_detection)->default_value(true),
			      "OUTLIER DETECTION FOR REAL ROBOT!")
		    ("rejection_multiplier", po::value<double>(&pflags.out_reject_mult)->default_value(2.0),
				"OUTLIER DETECTION MULTIPLIER FOR REAL ROBOT!")
			("check_bounce", po::value<bool>(&pflags.check_bounce),
			      "checking bounce before moving!")
		    ("weights", po::value<std::vector<double>>(&pflags.weights)->multitoken(), "hit,net,land weights for DP")
			("mult_vel", po::value<std::vector<double>>(&pflags.mult_vel)->multitoken(), "velocity mult. for DP")
			("penalty_loc", po::value<std::vector<double>>(&pflags.penalty_loc)->multitoken(), "punishment locations for DP")
			("rest_posture_optim", po::value<bool>(&pflags.optim_rest_posture)->default_value(false),
				"turn on resting state optimization")
			//("lookup", po::value<bool>(&flags.lookup), "start moving with lookup")
			("algorithm", po::value<int>(&alg_num)->default_value(0),
				  "optimization method")
			("mpc", po::value<bool>(&pflags.mpc)->default_value(false),
				 "corrections (MPC)")
			("spin", po::value<bool>(&pflags.spin)->default_value(false),
						 "apply spin model")
			("verbose", po::value<int>(&pflags.verbosity)->default_value(1),
		         "verbosity level")
		    ("save_data", po::value<bool>(&pflags.save)->default_value(false),
		         "saving robot/ball data")
		    ("ball_land_des_x_offset", po::value<double>(&pflags.ball_land_des_offset[0]),
		    	 "ball land x offset")
			("ball_land_des_y_offset", po::value<double>(&pflags.ball_land_des_offset[1]),
				 "ball land y offset")
		    ("time_land_des", po::value<double>(&pflags.time_land_des),
		    	 "time land des")
			("start_optim_offset", po::value<double>(&pflags.optim_offset),
				 "start optim offset")
			("time2return", po::value<double>(&pflags.time2return),
						 "time to return to start posture")
			("freq_mpc", po::value<int>(&pflags.freq_mpc), "frequency of updates")
		    ("min_obs", po::value<int>(&pflags.min_obs), "minimum obs to start filter")
		    ("var_noise", po::value<double>(&pflags.var_noise), "std of filter obs noise")
		    ("var_model", po::value<double>(&pflags.var_model), "std of filter process noise")
		    ("t_reset_threshold", po::value<double>(&pflags.t_reset_thresh), "filter reset threshold time")
		    ("VHPY", po::value<double>(&pflags.VHPY), "location of VHP")
		    ("url", po::value<std::string>(&pflags.zmq_url), "TCP URL for ZMQ connection")
		    ("debug_vision", po::value<bool>(&pflags.debug_vision), "print ball in listener")
			("listen_2d", po::value<bool>(&pflags.listen_2d), "listen to 2d server if true");
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
    set_player_algorithm(alg_num);
    pflags.detach = true; // always detached in SL/REAL ROBOT!
    return reinterpret_cast<void*>(&pflags);
}



void play(const SL_Jstate joint_state[NDOF+1],
          SL_DJstate joint_des_state[NDOF+1]) {

    // connect to ZMQ server
    // acquire ball info from ZMQ server
    // if new ball add status true else false
    // call play function
    static Listener* listener;

    // since old code support multiple blobs
    static ball_obs blob;
	static vec7 q0;
	static optim::joint qact;
	static optim::joint qdes;
	static Player *robot = nullptr;

	if (pflags.reset) {
		delete listener;
		listener = new Listener(pflags.zmq_url,
								pflags.listen_2d,
								pflags.debug_vision);
		for (int i = 0; i < NDOF; i++) {
			qdes.q(i) = q0(i) = joint_state[i+1].th;
			qdes.qd(i) = 0.0;
			qdes.qdd(i) = 0.0;
		}
		delete robot;
		robot = new Player(q0,pflags);
		pflags.reset = false;
	}
	else {
	    // update ball info
	    listener->fetch(blob);

		for (int i = 0; i < NDOF; i++) {
			qact.q(i) = joint_state[i+1].th;
			qact.qd(i) = joint_state[i+1].thd;
			qact.qdd(i) = joint_state[i+1].thdd;
		}
		robot->play(blob,qact,qdes);
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
    static Player *cp = nullptr; // centered player

    if (pflags.reset) {
        for (int i = 0; i < NDOF; i++) {
            qdes.q(i) = q0(i) = joint_state[i+1].th;
            qdes.qd(i) = 0.0;
            qdes.qdd(i) = 0.0;
        }
        cp = new Player(q0,pflags);
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

void save_joint_data(const SL_Jstate joint_state[NDOF+1],
          	  	  const SL_DJstate joint_des_state[NDOF+1],
				  const int save_qdes,
				  const int reset) {

	using namespace std::chrono;
	const std::string home = std::getenv("HOME");
    const std::string joint_file = home + "/table-tennis/joints.txt";

    static rowvec q = zeros<rowvec>(NDOF);
    static rowvec qdes = zeros<rowvec>(NDOF);
    static std::ofstream stream_joints;

    if (reset) {
    	stream_joints.close();
        stream_joints.open(joint_file,std::ofstream::out);
    }
    else {
    	milliseconds ms = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    	double time = ms.count();
		for (int i = 1; i <= NDOF; i++) {
			q(i-1) = joint_state[i].th;
			if (save_qdes) {
				qdes(i-1) = joint_des_state[i].th;
			}
		}

		if (stream_joints.is_open()) {
			if (save_qdes) {
				stream_joints << std::fixed << time << join_horiz(q,qdes);
			}
			else {
				stream_joints << std::fixed << time << q;
			}
		}
		//stream_joints.close();
    }
}

void save_ball_data(const char* url_string,
                     const int listen_2d,
                     const int debug_vision,
                     const int reset) {

	using namespace std::chrono;
    static Listener* listener;
    static ball_obs obs;
    static std::ofstream stream_balls;
    static std::string home = std::getenv("HOME");
    static std::string ball_file = home + "/table-tennis/balls.txt";

    if (reset) {
    	delete listener;
    	listener = new Listener(url_string,listen_2d,(bool)debug_vision);
    	stream_balls.close();
    	stream_balls.open(ball_file, std::ofstream::out);
    }
    else {
    	milliseconds ms = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    	double time = ms.count();
        // update blobs structure with new ball data
        listener->fetch(obs);
        if (obs.status) {
            if (stream_balls.is_open()) {
                stream_balls << std::fixed << time << obs.pos.t();
            }
        }
    }
}

void* load_serve_options(double custom_pose[], serve_task_options *options) {

    namespace po = boost::program_options;
    using std::string;
    const string home = std::getenv("HOME");
    string config_file = home + "/table-tennis/config/serve.cfg";

    try {
        // Declare a group of options that will be
        // allowed in config file
        po::options_description config("Configuration");
        config.add_options()
        ("use_rbf", po::value<bool>(&sflags.use_rbf),
        "Use Radial Basis Functions (RBF) instead of DMPs")
        ("json_file", po::value<string>(&sflags.json_file),
        "JSON File to load movement from")
        ("start_from_act_state", po::value<bool>(&sflags.start_from_act_state)->default_value(false),
        "Evolve from actual sensor state vs. desired state")
        ("save_act_joint", po::value<bool>(&sflags.save_joint_act_data),
        "Save act joint data during movement to txt file")
        ("save_des_joint", po::value<bool>(&sflags.save_joint_des_data),
        "Save des joint data during movement to txt file")
        ("save_ball_data", po::value<bool>(&sflags.save_ball_data),
        "Save ball data acquired via Listener to txt file")
        ("use_inv_dyn_fb", po::value<bool>(&sflags.use_inv_dyn_fb),
        "Use computed-torque control if false")
        ("url", po::value<std::string>(&sflags.zmq_url), "TCP URL for ZMQ connection")
        ("debug_vision", po::value<bool>(&sflags.debug_vision), "print ball in listener")
		("listen_2d", po::value<bool>(&sflags.listen_2d), "listen to 2d server if true")
        ("detach", po::value<bool>(&sflags.detach)->default_value(true),
              "detach optimization if true")
        ("mpc", po::value<bool>(&sflags.mpc)->default_value(false),
             "run optimization to correct for mispredictions etc.")
        ("freq_mpc", po::value<int>(&sflags.freq_mpc)->default_value(1),
             "frequency of running optimization for corrections")
        ("verbose", po::value<bool>(&sflags.verbose)->default_value(false),
             "verbosity level")
        ("ball_land_des_x_offset", po::value<double>(&sflags.ball_land_des_x_offset),
             "ball land x offset")
        ("ball_land_des_y_offset", po::value<double>(&sflags.ball_land_des_y_offset),
             "ball land y offset")
        ("time_land_des", po::value<double>(&sflags.time_land_des),
             "time land des");
        po::variables_map vm;
        std::ifstream ifs(config_file.c_str());
        if (!ifs) {
            cout << "can not open config file: " << config_file << "\n";
        }
        else {
            po::store(parse_config_file(ifs, config), vm);
            notify(vm);
        }
    }
    catch(std::exception& e) {
        cout << e.what() << "\n";
    }
    options->use_inv_dyn_fb = sflags.use_inv_dyn_fb;

    string json_file = home + "/table-tennis/json/" + sflags.json_file;
    vec7 pose;
    if (sflags.use_rbf) {
    	RBF rbfs = RBF(json_file);
    	rbfs.get_init_pos(pose);
    }
    else {
    	using dmps = Joint_DMPs;
        dmps multi_dmp = dmps(json_file);
        multi_dmp.get_init_pos(pose);
    }

    for (int i = 0; i < NDOF; i++) {
        custom_pose[i] = pose(i);
    }
    sflags.reset = true;
    return reinterpret_cast<void*>(&sflags);
}

void serve_ball(const SL_Jstate joint_state[],
                 SL_DJstate joint_des_state[]) {

    static Listener* listener;
    static ball_obs blob;
    static optim::joint qact;
    static optim::joint qdes;
    static ServeBall<dmps> *robot_starts_with_dmp = nullptr; // class for serving ball
    static ServeBall<RBF> *robot_starts_with_rbf = nullptr;

    if (sflags.reset) {
    	delete listener;
    	listener = new Listener(sflags.zmq_url,
    							sflags.listen_2d,
								sflags.debug_vision);
        vec7 qinit;
        if (sflags.start_from_act_state) {
            for (int i = 0; i < NDOF; i++) {
                qinit(i) = joint_state[i+1].th;
            }
        }
        for (int i = 0; i < NDOF; i++) {
            qdes.q(i) = joint_state[i+1].th;
            qdes.qd(i) = 0.0;
            qdes.qdd(i) = 0.0;
        }
        if (sflags.use_rbf) {
            delete robot_starts_with_rbf;
        	robot_starts_with_rbf = new ServeBall<RBF>(sflags,qinit);
        }
        else {
            delete robot_starts_with_dmp;
        	robot_starts_with_dmp = new ServeBall<dmps>(sflags,qinit);
        }
        sflags.reset = false;
    }
    else {
        // update blobs structure with new ball data
        listener->fetch(blob);
        for (int i = 0; i < NDOF; i++) {
            qact.q(i) = joint_state[i+1].th;
            qact.qd(i) = joint_state[i+1].thd;
            qact.qdd(i) = joint_state[i+1].thdd;
        }
        if (sflags.use_rbf) {
        	robot_starts_with_rbf->serve(blob,qact,qdes);
        }
        else {
        	robot_starts_with_dmp->serve(blob,qact,qdes);
        }
    }

    // update desired joint state
    for (int i = 0; i < NDOF; i++) {
        joint_des_state[i+1].th = qdes.q(i);
        joint_des_state[i+1].thd = qdes.qd(i);
        joint_des_state[i+1].thdd = qdes.qdd(i);
    }
}
