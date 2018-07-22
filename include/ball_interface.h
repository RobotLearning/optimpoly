#ifndef BALL_INTERF_H
#define BALL_INTERF_H

#include <vector>
#include <map>
#include "sl_structs.h"

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

#endif //BALL_INTERF_H
