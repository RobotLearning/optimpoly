#include <armadillo>
#include "constants.h"

using arma::zeros;
using arma::vec;
using arma::vec3;
using arma::vec7;
using arma::mat;

namespace DMP {

using namespace const_tt;

/**
 * @brief Desired/actual joint positions, velocities, accelerations.
 *
 * Main Player function play() updates desired joint values
 * every 2ms for Barrett WAM.
 */
struct joint {
    vec7 q = zeros<vec>(NDOF); //!< desired joint pos
    vec7 qd = zeros<vec>(NDOF); //!< desired joint vels
    vec7 qdd = zeros<vec>(NDOF); //!< desired joint accs
};

/**
 * \brief Canonical 1-D system in phase space to evolve (multiple) DMPs.
 */
struct Canonical {
    double a = 1.0;
    double tau = 1.0;
    double phase = 1.0;

    /** \brief Reset phase back to 1.0 */
    void reset();
    /** \brief Make a small step of dt sec. with Euler discretization */
    void step(const double & dt);
};

/**
 * @brief Discrete DMP for one degree of freedom
 */
class DMP {

private:
    vec3 x0 = zeros<vec>(3); //!< initial state of the DMP
    vec3 x = zeros<vec>(3); //!< state of the DMP
    vec w = zeros<vec>(10); //!< weights of the DMP
    vec h = zeros<vec>(10); //!< heights of the basis functions
    vec c = zeros<vec>(10); //!< centers of the basis functions
    double alpha = 25.0; //!< time constant
    double beta = 25/4.0; //!< time constant
    double g; //!< goal state

    /** \brief  Create a vector of basis functions evaluated at current phase x */
    vec basis(const double & x) const;

    /** \brief Create the (forces) accelerations due to the weights */
    double forcing(const double & phase) const;

public:

    /** \brief Initialize the DMP weights, goal, init. state etc. */
    DMP(std::vector<double> weights_,
        std::vector<double> centers_,
        std::vector<double> heights_,
        double alpha_,
        double beta_,
        double goal_,
        double init_pos);

    /** \brief Re-initialize the DMP */
    void reset();

    /** \brief Set weights, centers and heights of the basis functions */
    void set_weights(const std::vector<double> & w_stdvec,
                     const std::vector<double> & c_stdvec,
                     const std::vector<double> & h_stdvec);

    void set_goal_state(const double & goal);
    void get_goal_state(double & goal) const;
    void set_init_state(const double & init_pos);
    void get_init_state(double & init_pos) const;

    /** \brief Make a small step of dt seconds with Euler discretization */
    vec3 step(const Canonical & can, const double & dt);

};

/**
 * @brief Class holding one discrete DMP for each joint.
 *
 */
class Joint_DMPs {

private:
    Canonical can = Canonical(); //!< Canonical system with phase
    std::vector<DMP> dmps; //!< one discrete DMP for each joint

    /** @brief Reset the dmps and the canonical system */
    void reset();

public:

    /** \brief Empty constructor */
    Joint_DMPs() {};

    /** @brief Load DMPs from a json file. */
    Joint_DMPs(const std::string & param_file);

    /** @brief Evolve the DMPs dt seconds */
    void step(const double & dt, joint & Q);

    /** @brief Evolve the DMPs DT seconds for a total duration of T sec.*/
    mat evolve(const double & T);

    /** @brief Set initial positions of the DMPs */
    void set_init_pos(const vec & pos);

    /** @brief Get initial positions of the DMPs */
    void get_init_pos(vec & pos) const;

    /** @brief Set goal positions of the DMPs */
    void set_goal_pos(const vec & pos);

    /** @brief Return goal positions of the DMPs */
    void get_goal_pos(vec & pos) const;

    /** @brief Return speed of the movement */
    double get_time_constant() const;

    /** \brief Set speed of the movement */
    void set_time_constant(const double & tau);
};

}
