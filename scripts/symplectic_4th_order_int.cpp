void check_symplectic_int_4th() { /* It should be more accurate with fewer
                                     cycles */

  arma_rng::set_seed_random();
  double topspin = -50.0;
  TableTennis tt = TableTennis(true, false);
  tt.set_topspin(topspin);
  int ball_launch_side = (randi(1, distr_param(0, 2)).at(0));

  tt.set_ball_gun(0.05, ball_launch_side);
  vec6 init_state = tt.get_ball_state();
  int N = 25; // it should not bounce
  double dt = 0.02;
  for (int i = 0; i < N; i++) {
    tt.integrate_ball_state(dt);
  }
  vec6 ball_pred1 = tt.get_ball_state();

  tt.set_ball_state(init_state);
  for (int i = 0; i < 5; i++)
    tt.symplectic_int_fourth(0.1);
  // tt.symplectic_int_fourth(0.5);
  vec6 ball_pred2 = tt.get_ball_state();

  std::cout << ball_pred1.t();
  std::cout << ball_pred2.t();
}
