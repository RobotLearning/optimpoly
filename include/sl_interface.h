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

#include <map>
#include "constants.h"

/* The data structures from SL */
/** @brief (actual) joint space state for each DOF */
struct SL_Jstate {
    double   th;   /*!< theta */
    double   thd;  /*!< theta-dot */
    double   thdd; /*!< theta-dot-dot */
    double   ufb;  /*!< feedback portion of command */
    double   u;    /*!< torque command */
    double   load; /*!< sensed torque */
};
/** @brief (desired) joint space state commands for each DOF */
struct SL_DJstate { /*!< desired values for controller */
    double   th;   /*!< theta */
    double   thd;  /*!< theta-dot */
    double   thdd; /*!< theta-dot-dot */
    double   uff;  /*!< feedforward torque command */
    double   uex;  /*!< externally imposed torque */
};

/** @brief (actual) Cartesian state */
struct SL_Cstate {
    double x[const_tt::NCART+1];    /*!< Position [x,y,z] */
    double xd[const_tt::NCART+1];   /*!< Velocity */
    double xdd[const_tt::NCART+1];  /*!< Acceleration */
};

/** @brief Vision blob info coming from SL (after calibration). */
struct blob_state {
    int status = 0; //!< was ball detected reliably in cameras
    double pos[const_tt::NCART] = {0.0,0.0,0.0}; //!< ball center cartesian positions from cameras 1 and 2(after calibration)
};

/** @brief Listens continuously to ball info sent via ZMQ server. */
class Listener {

private:
    const std::string url; //!< TCP connection to computer and port
    std::map<double, std::vector<double>> obs; //!< time stamp as key and 3d vec as obs
    unsigned int max_obs_saved = 1000; //!< limit of observation map
    bool active = false; //!< actively listening the port
    bool debug = false; //!< for printing ZMQ connection info
    bool new_data = false; //!< new ball data has been saved but not fetched

    /** @brief Listen to TCP port broadcast via ZMQ based 3D ball server */
    void listen();

public:

    /** @brief Constructor that detaches the listen() private function in a thread. */
    Listener(const std::string & url = "tcp://helbe:7660", bool debug_ = false);

    /** @brief Stop listening. */
    void stop();

    /** @brief Fetch latest new ball data to blob structure */
    void fetch(blob_state & blob);

    /** @brief Returns number of observations saved, prints in debug mode*/
    int give_info();
};

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
