#include "ball_interface.h"
#include "json.hpp"
#include "zmqpp/zmqpp.hpp"
#include <thread>

using json = nlohmann::json;
using std::cout;
using std::endl;

void Listener::listen3d() {

    // initialize the zmq
    zmqpp::context context;
    zmqpp::socket_type type = zmqpp::socket_type::sub;
    zmqpp::socket socket(context,type);
    socket.set(zmqpp::socket_option::subscribe, "");
    socket.connect(url);

    if (debug)
        cout << "Starting listening to 3D ball pos..." << endl;

    while (active) {
        zmqpp::message msg;
        socket.receive(msg);
        std::string body;
        msg >> body;
        try {
            json jobs = json::parse(body);
            unsigned num = jobs.at("num");
            double time = jobs.at("time");
            json obs_j = jobs.at("obs");
            ball_pos obs_3d = {obs_j[0], obs_j[1], obs_j[2]};
            obs3d[num] = obs_3d;
            new_data = true;
            if (debug) {
                std::string ball = "[" +
                        std::to_string(obs_3d[0]) + " " +
                        std::to_string(obs_3d[1]) + " " +
                        std::to_string(obs_3d[2]) + "]";
                cout << "Received item " << ball << " at time: " << time << endl;
            }
            // keep size to a max
            if (obs3d.size() > max_obs_saved) {
                obs3d.erase(obs3d.begin());
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

void Listener::listen2d() {

    // initialize the zmq
    zmqpp::context context;
    zmqpp::socket_type type = zmqpp::socket_type::sub;
    zmqpp::socket socket(context,type);
    socket.set(zmqpp::socket_option::subscribe, "");
    socket.connect(url);

    if (debug)
        cout << "Starting listening to 2D pixel data..." << endl;

    while (active) {
        zmqpp::message msg;
        socket.receive(msg);
        std::string body;
        msg >> body;
        try {
            json jobs = json::parse(body);
            unsigned int num = jobs.at("num");
            double time = jobs.at("time");

            if (!jobs.at("obs").is_null()) {
            	pixels x;
            	x.cam_id = jobs.at("cam_id");
            	json jx = jobs.at("obs");
            	x.vals[0] = jx[0];
            	x.vals[1] = jx[1];
            	obs2d[num].push_back(x);
            	if (debug) {
            		std::string ball = "[" +
            				std::to_string(x.vals[0]) + " " +
							std::to_string(x.vals[1]) + "]";
            		cout << "Num: " << num << ". Received pixels" << ball << " for cam: " << x.cam_id << endl;
            	}
            }
        }
        catch (const std::exception & ex) {
            if (debug) {
                cout << "No pixels received..." << ex.what() << endl;
            }
        }
    }
    if (debug)
        cout << "Finished listening..." << endl;
}

void Listener::triangulate() {

	for (auto iter = obs2d.rbegin(); iter != obs2d.rend(); ++iter) {
		if (iter->second.size() >= 2) {
			unsigned num = iter->first;
			if (debug)
				cout << "Triangulating num: " << num<< endl;

			// triangulate by solving svd
			// TODO:
            ball_pos obs_3d;
            obs3d[num] = obs_3d;
            new_data = true;
            if (debug) {
                std::string ball = "[" +
                        std::to_string(obs_3d[0]) + " " +
                        std::to_string(obs_3d[1]) + " " +
                        std::to_string(obs_3d[2]) + "]";
                cout << "Received item " << ball << " at time: " << time << endl;
            }
            // keep size to a max
            if (obs3d.size() > max_obs_saved) {
                obs3d.erase(obs3d.begin());
            }
		}
	}
}

Listener::Listener(const std::string & url_,
					const bool run_2d,
					const bool debug_) : url(url_), debug(debug_) {
    using std::thread;
    active = true;
    if (run_2d) {
        thread t = thread(&Listener::listen2d,this);
        t.detach();
        thread t = thread(&Listener::triangulate,this);
        t.detach();
    }
    else {
        thread t = thread(&Listener::listen3d,this);
        t.detach();
    }

    //t.join();
}

void Listener::stop() {
    active = false;
    new_data = false;
}

void Listener::fetch(ball_obs & blob) { // fetch the latest data

    blob.status = false;
    // if there is a latest new data that was unread
    if (new_data) {
        blob.status = true;
        ball_pos latest_data = obs3d.rbegin()->second;
        blob.pos[0] = latest_data[0];
        blob.pos[1] = latest_data[1];
        blob.pos[2] = latest_data[2];
        new_data = false;
    }
}

int Listener::give_info() {

    using std::cout;

    if (debug) {
        cout << "Printing data...\n";
        cout << "==================================\n";
        cout << "Time \t obs[x] \t obs[y] \t obs[z]\n";
        for (std::pair<const double, std::vector<double>> element : obs3d) {
            cout << element.first << " \t"
                    << element.second[0] << " \t"
                    << element.second[1] << " \t"
                    << element.second[2] << "\n";
        }
    }
    return obs3d.size();
}
