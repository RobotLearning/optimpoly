#include <zmqpp/zmqpp.hpp>
#include <iostream>
#include <unordered_map>
#include "json.hpp"

using std::string;
using std::cout;
using std::endl;
using json = nlohmann::json;

int main() {
    const string url = "tcp://helbe:7660";

    // initialize the zmq
    zmqpp::context context;
    zmqpp::socket_type type = zmqpp::socket_type::sub;
    zmqpp::socket socket(context,type);
    socket.connect(url);
    socket.set(zmqpp::socket_option::subscribe, "");

    const int num_obs_max = 1000;
    std::unordered_map<double, std::vector<double>> obs;

    cout << "Starting listening..." << endl;
    while(obs.size() < num_obs_max) {
        zmqpp::message msg;
        socket.receive(msg);
        string body;
        msg >> body;
        json jobs = json::parse(body);
        unsigned int num = jobs.at("num");
        double time = jobs.at("time");
        json obs_j = jobs.at("obs");
        std::vector<double> obs_3d = {obs_j[0], obs_j[1], obs_j[2]};
        obs[time] = obs_3d;
        //cout << "Received item at time: " << time << endl;
    }
    cout << "Finished listening..." << endl;

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
