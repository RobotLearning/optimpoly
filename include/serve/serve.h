/*
 * serve.h
 *
 *  Created on: Jul 22, 2018
 *      Author: okoc
 */

#ifndef INCLUDE_SERVE_H_
#define INCLUDE_SERVE_H_

namespace serve {

using dmps = Joint_DMPs;
using vec_str = std::vector<std::string>;

struct serve_flags {
    bool reset = false; //!< reset the serve class
    bool save_joint_act_data = false;
    bool save_joint_des_data = false;
    bool save_ball_data = false;
    bool start_dmp_from_act_state = false;
    bool use_inv_dyn_fb = false; //!< in SL apply inv. dyn. feedback
    std::string json_file = "dmp4.json"; //!< json file to load dmp from
    std::string zmq_url = "tcp://helbe:7660"; //!< URL for ZMQ connection
    bool debug_vision = false; //!< print received vision info
};

class ServeBall {

private:
    serve_flags sflags;
    double T = 1.0;
    vec7 q_rest_des;
    vec7 q_hit_des;
    optim::Optim *opt = nullptr; // optimizer

public:

    /**
     * \brief Serve a table tennis ball
     *
     * Initialize the serve with a DMP and then repeatedly correct
     * with optimization whenever the predicted ball is not hit by the
     * predicted movement.
     */
    void serve(const player::EKF & filter,
                dmps & multi_dmp,
                joint & Q);

    ServeBall(const double & T, dmps & multi_dmp);
    ~ServeBall();

};

/** \brief Initialize DMP from a random JSON file */
dmps init_dmps();

/** \brief Estimate ball state using a few initial ball estimates.*/
void estimate_ball_state(const vec3 & obs, player::EKF & filter);

/** \brief Utility function to return vector of JSON files from subfolder */
vec_str get_files(std::string folder_name);

}

#endif /* INCLUDE_SERVE_SERVE_H_ */
