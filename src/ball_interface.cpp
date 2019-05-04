#include "ball_interface.h"
#include "json.hpp"
#include "utils.h"
#include "zmqpp/zmqpp.hpp"
#include <thread>

using json = nlohmann::json;
using std::endl;

void Listener::listen3d() {

  // initialize the zmq
  zmqpp::context context;
  zmqpp::socket_type type = zmqpp::socket_type::sub;
  zmqpp::socket socket(context, type);
  socket.set(zmqpp::socket_option::subscribe, "");
  socket.connect(url_);

  if (debug_) {
    stream_balls_ << "Starting listening to 3D ball pos..." << endl;
  }

  while (active_) {
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
      obs3d_[num] = obs_3d;
      new_data_ = true;
      if (debug_) {
        std::string ball = "[" + std::to_string(obs_3d[0]) + " " +
                           std::to_string(obs_3d[1]) + " " +
                           std::to_string(obs_3d[2]) + "]";
        stream_balls_ << "Received item " << ball << " at time: " << time
                      << endl;
      }
      // keep size to a max
      if (obs3d_.size() > max_obs_saved_) {
        obs3d_.erase(obs3d_.begin());
      }
    } catch (const std::exception &ex) {
      if (debug_) {
        stream_balls_ << "No ball detected..." << ex.what() << endl;
      }
    }
  }
  if (debug_) {
    stream_balls_ << "Finished listening..." << endl;
  }
}

void Listener::listen2d() {

  // initialize the zmq
  zmqpp::context context;
  zmqpp::socket_type type = zmqpp::socket_type::sub;
  zmqpp::socket socket(context, type);
  socket.set(zmqpp::socket_option::subscribe, "");
  socket.connect(url_);

  if (debug_) {
    stream_balls_ << "Starting listening to 2D pixel data..." << endl;
  }

  while (active_) {
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
        x.time_stamp = time;
        x.cam_id = jobs.at("cam_id");
        json jx = jobs.at("obs");
        x.vals[0] = jx[0];
        x.vals[1] = jx[1];
        obs2d_[num].push_back(x);
        if (debug_) {
          std::string ball = "[" + std::to_string(x.vals[0]) + " " +
                             std::to_string(x.vals[1]) + "]";
          stream_balls_ << "Num: " << num << ". Received pixels" << ball
                        << " for cam: " << x.cam_id << endl;
        }
      }
      // keep size to a max
      if (obs2d_.size() > max_obs_saved_) {
        obs2d_.erase(obs2d_.begin());
      }
    } catch (const std::exception &ex) {
      if (debug_) {
        stream_balls_ << "No pixels received..." << ex.what() << endl;
      }
    }
  }
  if (debug_) {
    stream_balls_ << "Finished listening..." << endl;
  }
}

void Listener::convert_to_3d() {

  if (debug_) {
    stream_balls_ << "Starting triangulating from 2D pixels to 3D pos..."
                  << endl;
  }

  while (active_) {
    for (auto iter = obs2d_.cbegin(); iter != obs2d_.cend();
         /* no increment */) {
      if (iter->second.size() >= 2) {
        unsigned num = iter->first;
        if (debug_) {
          stream_balls_ << "Triangulating num: " << num << endl;
        }

        // triangulate by solving svd
        ball_pos obs_3d;
        bool success =
            triangulate(calib_mats_, iter->second, triangulation_, obs_3d);
        if (success) {
          obs2d_.erase(iter++);
          obs3d_[num] = obs_3d;
          new_data_ = true;
        }
        if (debug_) {
          std::string ball = "[" + std::to_string(obs_3d[0]) + " " +
                             std::to_string(obs_3d[1]) + " " +
                             std::to_string(obs_3d[2]) + "]";
          stream_balls_ << "Received item " << ball << " at num: " << num
                        << endl;
        }
        // keep size to a max
        if (obs3d_.size() > max_obs_saved_) {
          obs3d_.erase(obs3d_.begin());
        }
      } else {
        ++iter;
      }
    }
  }
}

Listener::Listener(const std::string &url, const bool run_2d, const bool debug,
                   const std::string triangulation)
    : url_(url), debug_(debug), triangulation_(triangulation),
      max_obs_saved_(1000), active_(true), new_data_(false) {
  using std::thread;
  std::string home = std::getenv("HOME");
  std::string debug_file = home + "/projects/table-tennis/debug_listener.txt";
  stream_balls_.open(debug_file, std::ofstream::out);
  active_ = true;
  if (run_2d) {
    // stream_balls_ << "Launching 2D listener...\n";
    calib_mats_ = load_proj_mats("server_3d_conf_ping_okan.json");
    thread t = thread(&Listener::listen2d, this);
    t.detach();
    thread t2 = thread(&Listener::convert_to_3d, this);
    t2.detach();
  } else {
    // stream_balls_ << "Launching 3D listener...\n";
    thread t = thread(&Listener::listen3d, this);
    t.detach();
  }
}

Listener::~Listener() { stream_balls_.close(); }

void Listener::stop() {
  active_ = false;
  new_data_ = false;
  stream_balls_.close();
  obs2d_.clear();
  obs3d_.clear();
}

void Listener::fetch(ball_obs &blob) { // fetch the latest data

  blob.status = false;
  // if there is a latest new data that was unread
  if (new_data_) {
    blob.status = true;
    blob.pos = obs3d_.rbegin()->second;
    new_data_ = false;
  }
}

int Listener::give_info() {

  if (debug_) {
    stream_balls_ << "Printing data...\n";
    stream_balls_ << "==================================\n";
    stream_balls_ << "Time \t obs[x] \t obs[y] \t obs[z]\n";
    for (std::pair<const double, ball_pos> element : obs3d_) {
      stream_balls_ << element.first << " \t" << element.second.t();
    }
  }
  return obs3d_.size();
}

/**
 * @brief Load projection matrices from json file
 */
std::map<unsigned, mat34>
load_proj_mats(const std::string &json_file = "server_3d_conf_ping.json") {

  using namespace arma;
  std::map<unsigned, mat34> calib_mats_;
  std::string home = std::getenv("HOME");
  std::string file_path = home + "/projects/table-tennis/json/" + json_file;
  std::ifstream stream(file_path);
  json jobs;
  stream >> jobs;
  json jcalib = jobs["stereo"]["calib"];
  for (auto elem : jcalib) {
    unsigned int id = elem.at("ID");
    mat34 val = json2mat(elem.at("val"));
    calib_mats_[id] = val;
  }
  return calib_mats_;
}

/**
 * Triangulate cameras 0 and 1, or cameras 2 and 3
 */
bool triangulate(const std::map<unsigned, mat34> &calib_mats_,
                 const std::vector<pixels> &obs_2d,
                 const std::string triangulation_method, ball_pos &pos3d) {

  const int NUM_CAMS = 4; // calib_mats_.size();
  const int NUM_PAIRS = NUM_CAMS / 2;
  using namespace arma;
  double pixels[NUM_CAMS][2];
  bool found[NUM_CAMS] = {false};

  for (auto p : obs_2d) {
    unsigned int idx = p.cam_id;
    pixels[idx][0] = p.vals[0];
    pixels[idx][1] = p.vals[1];
    found[idx] = true;
  }

  for (int i = 0; i < NUM_PAIRS; i++) {
    if (found[2 * i] && found[2 * i + 1]) {
      mat P1 = calib_mats_.at(2 * i);
      mat P2 = calib_mats_.at(2 * i + 1);

      // make sure mats are normalized
      /*P1 /= P1(2,3);
      P2 /= P2(2,3);*/

      vec4 pos4d;
      // METHOD ONE: Ax = 0
      if (triangulation_method.compare("DLT") == 0) {
        mat44 A = zeros<mat>(4, 4);
        A.row(0) = pixels[2 * i][0] * P1.row(2) - P1.row(0);
        A.row(1) = pixels[2 * i][1] * P1.row(2) - P1.row(1);
        A.row(2) = pixels[2 * i + 1][0] * P2.row(2) - P2.row(0);
        A.row(3) = pixels[2 * i + 1][1] * P2.row(2) - P2.row(1);

        mat44 U = zeros<mat>(4, 4);
        vec4 s = zeros<vec>(4);
        mat44 V = zeros<mat>(4, 4);
        svd(U, s, V, A);
        pos4d = V.tail_cols(1);
      } else { // METHOD TWO: invert P1 and P2
        mat44 P = join_vert(P1.rows(0, 1), P2.rows(0, 1));
        vec4 v_pixels = zeros<vec>(4);
        v_pixels(0) = pixels[2 * i][0];
        v_pixels(1) = pixels[2 * i][1];
        v_pixels(2) = pixels[2 * i + 1][0];
        v_pixels(3) = pixels[2 * i + 1][1];
        pos4d = solve(P, v_pixels);
      }
      pos3d = pos4d.head(3);
      pos3d /= pos4d(3); // normalize
      return true;
    }
  }
  return false;
}
