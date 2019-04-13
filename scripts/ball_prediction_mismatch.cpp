/*
 * Test mismatch in ball prediction
 * by using a spin model and a normal drag model for prediction
 *
 */
void check_mismatch_pred() {

	BOOST_TEST_MESSAGE("Checking ball prediction accuracy in mismatch case (spin)...");
	TableTennis tt = TableTennis(true, false); //init with spin
	arma_rng::set_seed(5);
	EKF filter = init_ball_filter(0.3,0.001);
	int N = 100;
	tt.set_ball_gun(0.2);
	vec6 ball_init = tt.get_ball_state();
	filter.set_prior(ball_init,eye<mat>(6,6));
	vec pos_errs, vel_errs;
	pos_errs = vel_errs = zeros<vec>(N);
	for (int i = 0; i < N; i++) {
		tt.integrate_ball_state(DT);
		filter.predict(DT,true);
		filter.update(tt.get_ball_position());
		pos_errs(i) = norm(filter.get_mean()(span(X,Z)) - tt.get_ball_position());
		vel_errs(i) = norm(filter.get_mean()(span(DX,DZ)) - tt.get_ball_velocity());
	}
	BOOST_TEST_MESSAGE("Pos error max: " << max(pos_errs));
	BOOST_TEST_MESSAGE("Vel error max: " << max(vel_errs));

}
