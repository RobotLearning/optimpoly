#include "ball_interface.h"
#include "json.hpp"
#include "zmqpp/zmqpp.hpp"
#include <thread>

using json = nlohmann::json;
using std::cout;
using std::endl;

static arma::mat json2mat(const json & jobj);
static std::map<unsigned, mat34> load_calib_mats(const std::string & json_file);
static ball_pos triangulate(const std::map<unsigned, mat34> & calib_mats,
							const std::vector<pixels> & obs_2d);

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
            // keep size to a max
            if (obs2d.size() > max_obs_saved) {
                obs2d.erase(obs2d.begin());
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

void Listener::convert_to_3d() {

	for (auto iter = obs2d.rbegin(); iter != obs2d.rend(); ++iter) {
		if (iter->second.size() >= 2) {
			unsigned num = iter->first;
			if (debug)
				cout << "Triangulating num: " << num<< endl;

			// triangulate by solving svd
            ball_pos obs_3d = triangulate(calib_mats, iter->second);
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
    	calib_mats = load_calib_mats("json/server_3d_conf_ping_okan.json");
        thread t = thread(&Listener::listen2d,this);
        t.detach();
        thread t2 = thread(&Listener::convert_to_3d,this);
        t2.detach();
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
        blob.pos = obs3d.rbegin()->second;
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
                 << element.second << endl;
        }
    }
    return obs3d.size();
}

/**
 * @brief Load calibration matrices from json file
 */
static std::map<unsigned, mat34> load_calib_mats(const std::string & json_file =
												"json/server_3d_conf_ping.json") {

	using namespace arma;
	std::map<unsigned,mat34> calib_mats;
    std::ifstream stream(json_file);
    json jobs;
    stream >> jobs;
    json jcalib = jobs["stereo"]["calib"];
    for (auto elem : jcalib) {
        unsigned int id = elem.at("ID");
        mat34 val = json2mat( elem.at("val") );
        calib_mats[id] = val;
    }
    return calib_mats;
}

static arma::mat json2mat(const json & jobj) {


	int n = jobj.size();
	int m = jobj[0].size();
	arma::mat arma_mat(n,m,arma::fill::zeros);

	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++) {
			arma_mat(i,j) = jobj[i][j];
		}
	}
	return arma_mat;
}

/**
 * Triangulate cameras 0 and 1
 * TODO:
 */
static ball_pos triangulate(const std::map<unsigned, mat34> & calib_mats,
							const std::vector<pixels> & obs_2d) {

	using namespace arma;
	ball_pos pos3d;
	mat34 P0, P1;
	bool zero_found, one_found = false;
	double pixels0[2], pixels1[2];
	for (auto p : obs_2d) {
		if (p.cam_id == 0) {
			P0 = calib_mats[0];
			pixels0[0] = p.vals[0];
			pixels0[1] = p.vals[1];
			zero_found = true;
		}
		else if (p.cam_id == 1) {
			P1 = calib_mats[1];
			pixels1[0] = p.vals[0];
			pixels1[1] = p.vals[1];
			one_found = true;
		}
	}
	if (zero_found && one_found) {
		mat44 A;
		A.row(0) = pixels0[0]*P0.row(2) - P0.row(0);
		A.row(1) = pixels0[1]*P0.row(2) - P0.row(1);
		A.row(2) = pixels1[0]*P1.row(2) - P1.row(0);
		A.row(3) = pixels1[1]*P1.row(2) - P1.row(1);
		mat U;
		vec s;
		mat V;
		svd(U,s,V,A);
		vec4 pos4d = V.tail_cols(0);
		pos3d = pos4d.head(3);
		pos3d /= pos4d(3);
	}
	return pos3d;
}
