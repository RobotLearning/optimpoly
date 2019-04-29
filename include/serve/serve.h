#pragma once
#include "ball_interface.h"
#include "dmp.h"

namespace serve {

using dmps = Joint_DMPs;
using vec_str = std::vector<std::string>;

/**
 * \brief Flags used in the SERVE class and also in the SERVE task in SL.
 */
struct serve_flags {
  bool detach = false;  //!< detach optimization
  bool mpc = true;      //!< run optimization if MP is predicted to miss target
  bool verbose = false; //!< print optim info if true
  bool reset = false;   //!< reset the serve class
  bool save_joint_act_data = false;
  bool save_joint_des_data = false;
  bool save_ball_data = false;
  bool start_from_act_state = false;
  bool use_inv_dyn_fb = false; //!< in SL apply inv. dyn. feedback
  bool debug_vision = false;   //!< print received vision info
  int alg = 0;                 //!< 0 = DMP, 1 = RBF, 2 = CRBF
  int freq_mpc = 1.0;          //!< how many times per minute to re-run optim
  double time_land_des = 0.6;  //!< desired time to land on robot court first
    std::vector<double> ball_land_des_offset = {0.0, 0.0}; //!<
  double speedup = 1;                  //!< for the DMP speed up
  bool listen_2d = true;               //!< listen to 2d server or 3d server
  std::string json_file = "dmp4.json"; //!< json file to load dmp from
  std::string zmq_url = "tcp://localhost:7650"; //!< URL for ZMQ connection
  std::string triangulation = "DLT";
};

template <typename T> class ServeBall {

private:
  bool ran_optim = false;
  double t_clock = 0.0;
  serve_flags sflags;
  arma::wall_clock timer;
  optim::spline_params p;
  T mp = T(); //!< movement pattern
  double Tmax = 1.0;
  vec7 q_rest_des;
  vec7 q_hit_des;
  int idx_balls_obs_filter = 0;
  bool init_filter = false;
  optim::Optim *opt = nullptr; // optimizer
  player::EKF filter = player::init_ball_filter(0.3, 0.001, false);

  /*
   * \brief Correct movement by running an optimization.
   * The plan is to serve the ball to a desired position on robot court.
   */
  void correct_with_optim(const joint &qact);

  /*
   * \brief If the filter was initialized then we have enough info to predict
   * what will happen to the ball. If a hit is predicted then returns true.
   */
  bool predict_ball_hit(const double &t_pred);

  /** \brief Estimate ball state using a few initial ball estimates.*/
  void estimate_ball_state(const ball_obs &obs);

public:
  /**
   * \brief Serve a table tennis ball
   *
   * Start by following a DMP and then repeatedly correct
   * with optimization whenever the predicted ball is not hit by the
   * predicted movement.
   */
  void serve(const ball_obs &obs, const joint &qact, joint &qdes);

  /**\brief Initialize serve class with a learned movement and setup
   * optimization for corrections later. */
  ServeBall(const serve_flags &sflags, const vec7 &qinit = zeros<vec>(7));
  ~ServeBall();

  void set_flags(const serve_flags &sflags);
};

/** \brief Initialize DMP from a random JSON file */
dmps init_dmps(std::string &file);

/** \brief Utility function to return vector of JSON files from subfolder */
vec_str get_files(std::string folder_name, std::string prefix);

} // namespace serve
