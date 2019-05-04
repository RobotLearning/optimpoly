#pragma once
#include <armadillo>
#include <map>
#include <vector>

using ball_pos = arma::vec3;
using mat34 = arma::mat::fixed<3, 4>;

struct pixels {
  double time_stamp = 0.0;    //!< time stamp of image
  unsigned cam_id = 0;        //!< camera ID
  std::array<double, 2> vals; //!< pixel values
};

/** @brief Vision ball info coming from vision computer (after calibration). */
struct ball_obs {
  bool status =
      false;      //!< was ball detected reliably in cameras and is it new data?
  arma::vec3 pos; //!< ball center cartesian positions from two cameras
};

/** @brief Listens continuously to ball info sent via ZMQ server. */
class Listener {

private:
   //!< stream for outputting debug info
  std::ofstream stream_balls_;
  std::string url_;            //!< TCP connection to computer and port
  std::map<unsigned, mat34>
      calib_mats_; //!< Calibration matrices loaded from json file
  std::map<unsigned, std::vector<pixels>, std::greater<unsigned>>
      obs2d_;                          //!< map with frame as key and pixels as
  std::map<unsigned, ball_pos> obs3d_; //!< frame id as key and 3d vec as obs

  unsigned const int max_obs_saved_ = 1000; //!< limit of observation map
  bool active_ = false;               //!< actively listening the port
  bool debug_ = false;                //!< for printing ZMQ connection info
  bool new_data_ = false; //!< new ball data has been saved but not fetched
  std::string triangulation_; //!< triangulation method, DLT or invert
  /** @brief Listen to TCP port broadcast via ZMQ based 3D ball server */
  void listen3d();

  /** @brief Listen to TCP port broadcast via ZMQ based 2D ball server */
  void listen2d();

  /** @brief Triangulate to 3d if at least 2 cameras have reported for a given
   * frame */
  void convert_to_3d();

public:
  /** @brief Constructor that detaches the listen() private function in a
   * thread.*/
  Listener(const std::string &url, bool run_2d, bool debug_, const std::string triangulation);
  ~Listener();

  /** @brief Stop listening. */
  void stop();

  /** @brief Fetch latest new ball data to blob structure */
  void fetch(ball_obs &obs);

  /** @brief Returns number of observations saved, prints in debug mode*/
  int give_info();
};

std::map<unsigned, mat34> load_proj_mats(const std::string &json_file);

/** @brief Triangulate from two 2d pixels to one 3d ball position. */
bool triangulate(const std::map<unsigned, mat34> &calib_mats,
                 const std::vector<pixels> &obs_2d,
                 const std::string triangulation,
                 ball_pos &obs_3d);
