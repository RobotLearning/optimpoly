void count_land_mpc() {

  BOOST_TEST_MESSAGE("Running MPC Test...");
  BOOST_TEST_MESSAGE("Counting ball landing with 3 different optim...");
  algo algs[] = {FOCUS, DP, VHP};

  for (int i = 0; i < 3; i++) {
    double Tmax = 1.0, lb[2 * NDOF + 1], ub[2 * NDOF + 1];
    set_bounds(lb, ub, 0.01, Tmax);
    vec7 lbvec(lb);
    vec7 ubvec(ub);
    TableTennis tt = TableTennis(true, true);
    int num_trials = 1;
    int num_lands = 0;
    int num_misses = 0;
    int num_not_valid = 0;
    arma_rng::set_seed_random();
    // arma_rng::set_seed(0);
    double std_noise = 0.0001;
    joint qact;
    vec3 obs;
    player_flags flags;
    flags.alg = algs[i];
    flags.mpc = true;
    flags.freq_mpc = 1;
    flags.verbosity = 0;
    Player *robot;
    int N = 2000;
    joint qdes = qact;
    racket robot_racket;
    int ball_launch_side;
    int joint_init_pose;

    for (int n = 0; n < num_trials; n++) { // for each trial
      std::cout << "Trial: " << n + 1 << std::endl;
      ball_launch_side = (randi(1, distr_param(0, 2)).at(0));
      joint_init_pose = (randi(1, distr_param(0, 2)).at(0));
      init_posture(qact.q, joint_init_pose, true);
      robot = new Player(qact.q, flags);
      tt.reset_stats();
      tt.set_ball_gun(0.05, ball_launch_side);
      // robot.reset_filter(std_model,std_noise);
      for (int i = 0; i < N; i++) { // one trial
        obs = tt.get_ball_position() + std_noise * randn<vec>(3);
        ball_obs obs_str;
        obs_str.status = true;
        obs_str.pos = obs;
        robot->play(obs_str, qact, qdes);
        // robot->cheat(qact, obs, qdes);
        calc_racket_state(qdes, robot_racket);
        tt.integrate_ball_state(robot_racket, DT);
        qact.q = qdes.q;
        qact.qd = qdes.qd;
      }
      if (tt.has_legally_landed()) {
        num_lands++;
      } else if (!tt.has_legally_bounced())
        num_not_valid++;
      else
        num_misses++;
      delete robot;
    }

    std::cout << "======================================================"
              << endl;
    std::cout << "Out of " << num_trials << " trials, " << num_lands
              << " lands, " << num_not_valid << "not valid balls, "
              << num_misses << " misses!" << std::endl;
    std::cout << "======================================================"
              << endl;
  }
}

void count_land() {

  BOOST_TEST_MESSAGE("Counting Robot Ball Landing with 3 different optim...");
  algo algs[] = {FOCUS, DP, VHP};
  for (int n = 0; n < 3; n++) {
    double Tmax = 2.0;
    double lb[2 * NDOF + 1], ub[2 * NDOF + 1];
    set_bounds(lb, ub, 0.01, Tmax);
    vec7 lbvec(lb);
    vec7 ubvec(ub);
    TableTennis tt = TableTennis(false, true);
    arma_rng::set_seed_random();
    // arma_rng::set_seed(5);
    tt.set_ball_gun(0.05, 0); // init ball on the centre
    double std_obs = 0.0001;  // std of the noisy observations
    joint qact;
    init_posture(qact.q, 1, false);
    vec3 obs;
    player_flags flags;
    flags.verbosity = 0;
    flags.mpc = false;
    flags.alg = algs[n];
    Player robot = Player(qact.q, flags);
    int N = 2000;
    joint qdes;
    qdes.q = qact.q;
    racket robot_racket;
    mat Qdes = zeros<mat>(NDOF, N);
    for (int i = 0; i < N; i++) {
      obs = tt.get_ball_position() + std_obs * randn<vec>(3);
      ball_obs obs_str;
      obs_str.status = true;
      obs_str.pos = obs;
      robot.play(obs_str, qact, qdes);
      // robot.cheat(qact, tt.get_ball_state(), qdes);
      Qdes.col(i) = qdes.q;
      calc_racket_state(qdes, robot_racket);
      // cout << "robot ball dist\t" << norm(robot_racket.pos -
      // tt.get_ball_position()) << endl;
      tt.integrate_ball_state(robot_racket, DT);
      // usleep(DT*1e6);
      qact.q = qdes.q;
      qact.qd = qdes.qd;
    }
    // cout << max(Qdes,1) << endl;

    vec2 ball_des;
    double time_des;
    robot.get_strategy(ball_des, time_des);
    BOOST_TEST_MESSAGE("Desired ball land: " << ball_des.t());
    BOOST_TEST_MESSAGE("Testing joint limits as well...");
    BOOST_TEST(all(max(Qdes, 1) < ubvec));
    BOOST_TEST(all(min(Qdes, 1) > lbvec));
    BOOST_TEST(tt.has_legally_landed());
    std::cout << "******************************************************"
              << std::endl;
  }
}

static void init_posture(vec7 &q0, int posture, bool verbose) {

  rowvec qinit;
  switch (posture) {
  case 2: // right
    if (verbose)
      cout << "Initializing robot on the right side.\n";
    qinit << 1.0 << -0.2 << -0.1 << 1.8 << -1.57 << 0.1 << 0.3 << endr;
    break;
  case 1: // center
    if (verbose)
      cout << "Initializing robot on the center.\n";
    qinit << 0.0 << 0.0 << 0.0 << 1.5 << -1.75 << 0.0 << 0.0 << endr;
    break;
  case 0: // left
    if (verbose)
      cout << "Initializing robot on the left side\n";
    qinit << -1.0 << 0.0 << 0.0 << 1.5 << -1.57 << 0.1 << 0.3 << endr;
    break;
  default: // default is the right side
    qinit << 1.0 << -0.2 << -0.1 << 1.8 << -1.57 << 0.1 << 0.3 << endr;
    break;
  }
  q0 = qinit.t();
}
