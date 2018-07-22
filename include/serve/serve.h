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

class ServeBall {

private:
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
