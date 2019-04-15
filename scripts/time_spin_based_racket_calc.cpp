void check_speed_spin_based_racket_calc() {
  // Testing speed of racket calculations with spin-based BVP solver
  // (optimization)

  arma_rng::set_seed_random();
  double topspin = -50.0;
  TableTennis tt = TableTennis(true, false);
  tt.set_topspin(topspin);
  int ball_launch_side = (randi(1, distr_param(0, 2)).at(0));
  tt.set_ball_gun(0.05, ball_launch_side);
  int N = 50;
  double dt = 0.02;
  mat balls_pred = zeros<mat>(6, N);
  for (int i = 0; i < N; i++) {
    tt.integrate_ball_state(dt);
    balls_pred.col(i) = tt.get_ball_state();
  }
  vec3 ball_land_des = {0.0, dist_to_table - 3 * table_length / 4.0,
                        floor_level - table_height + ball_radius};
  double time_land_des = 0.6;

  optim_des racket_des;
  racket_des.Nmax = N;
  static wall_clock timer;
  timer.tic();
  calc_racket_strategy(balls_pred, ball_land_des.head(2), time_land_des,
                       racket_des);
  double time1 = timer.toc() * 1e3;
  std::cout << "Elapsed time in ms: " << time1 << std::endl;
  timer.tic();
  calc_spin_racket_strategy(balls_pred, topspin, ball_land_des, time_land_des,
                            racket_des);
  double time2 = timer.toc() * 1e3;
  std::cout << "Elapsed time in ms: " << time2 << std::endl;
  // BOOST_TEST(time2 < time1 * 100);
}
