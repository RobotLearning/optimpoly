#ifndef BALL_INTERF_H
#define BALL_INTERF_H

#include <vector>
#include <map>
#include <armadillo>

using ball_pos = std::array<double,3>;

struct pixels {
	double time_stamp = 0.0; //!< time stamp of image
	unsigned cam_id = 0; //!< camera ID
	std::array<double,2> vals; //!< pixel values
};

/** @brief Vision ball info coming from vision computer (after calibration). */
struct ball_obs {
    bool status = false; //!< was ball detected reliably in cameras and is it new data?
    arma::vec3 pos; //!< ball center cartesian positions from two cameras
};

/** @brief Listens continuously to ball info sent via ZMQ server. */
class Listener {

private:
    const std::string url; //!< TCP connection to computer and port

    std::map<unsigned, std::vector<pixels>> obs2d; //!< map with frame as key and pixels as
    std::map<unsigned, ball_pos> obs3d; //!< frame id as key and 3d vec as obs

    unsigned int max_obs_saved = 1000; //!< limit of observation map
    bool active = false; //!< actively listening the port
    bool debug = false; //!< for printing ZMQ connection info
    bool new_data = false; //!< new ball data has been saved but not fetched

    /** @brief Listen to TCP port broadcast via ZMQ based 3D ball server */
    void listen3d();

    /** @brief Listen to TCP port broadcast via ZMQ based 2D ball server */
    void listen2d();

    /** @brief Triangulate to 3d if at least 2 cameras have reported for a given frame */
    void triangulate();

public:

    /** @brief Constructor that detaches the listen() private function in a thread.
     * TCP port is 7660 Helbe (vision computer) for 3d server
     *             7650 Helbe for 2d server
     *  */
    Listener(const std::string & url = "tcp://helbe:7660", bool run_2d, bool debug_ = false);

    /** @brief Stop listening. */
    void stop();

    /** @brief Fetch latest new ball data to blob structure */
    void fetch(ball_obs & obs);

    /** @brief Returns number of observations saved, prints in debug mode*/
    int give_info();
};

#endif //BALL_INTERF_H
