#include "ball_interface.h"
#include "json.hpp"
#include "zmqpp/zmqpp.hpp"

void Listener::listen() {

    using json = nlohmann::json;
    using std::cout;
    using std::endl;

    // initialize the zmq
    zmqpp::context context;
    zmqpp::socket_type type = zmqpp::socket_type::sub;
    zmqpp::socket socket(context,type);
    socket.set(zmqpp::socket_option::subscribe, "");
    socket.connect(url);

    if (debug)
        cout << "Starting listening..." << endl;

    while (active) {
        zmqpp::message msg;
        socket.receive(msg);
        std::string body;
        msg >> body;
        try {
            json jobs = json::parse(body);
            //unsigned int num = jobs.at("num");
            double time = jobs.at("time");
            json obs_j = jobs.at("obs");
            std::vector<double> obs_3d = {obs_j[0], obs_j[1], obs_j[2]};
            obs[time] = obs_3d;
            new_data = true;
            if (debug) {
                std::string ball = "[" +
                        std::to_string(obs_3d[0]) + " " +
                        std::to_string(obs_3d[1]) + " " +
                        std::to_string(obs_3d[2]) + "]";
                cout << "Received item " << ball << " at time: " << time << endl;
            }
            // keep size to a max
            if (obs.size() > max_obs_saved) {
                obs.erase(obs.begin());
            }
        }
        catch (const std::exception & ex) {
            if (debug) {
                cout << "No ball detected..." << ex.what() << endl;
            }
        }
    }
    if (debug)
        cout << "Finished listening..." << endl;
}

Listener::Listener(const std::string & url_, const bool debug_) : url(url_), debug(debug_) {
    using std::thread;
    active = true;
    thread t = thread(&Listener::listen,this);
    t.detach();
    //t.join();
}

void Listener::stop() {
    active = false;
    new_data = false;
}

void Listener::fetch(blob_state & blob) { // fetch the latest data

    blob.status = false;
    // if there is a latest data
    if (new_data) {
        blob.status = true;
        std::vector<double> last_data = obs.rbegin()->second;
        blob.pos[0] = last_data[0];
        blob.pos[1] = last_data[1];
        blob.pos[2] = last_data[2];
        new_data = false;
    }
}

int Listener::give_info() {

    using std::cout;

    if (debug) {
        cout << "Printing data...\n";
        cout << "==================================\n";
        cout << "Time \t obs[x] \t obs[y] \t obs[z]\n";
        for (std::pair<const double, std::vector<double>> element : obs) {
            cout << element.first << " \t"
                    << element.second[0] << " \t"
                    << element.second[1] << " \t"
                    << element.second[2] << "\n";
        }
    }
    return obs.size();
}
