/**
 * Evaluate time efficiency of FP and DP optimizations over hundreds of tests
 */
void test_time_efficiency() {

	// For FP or DP
	// initialize ball randomly each time
	// measure time elapsed for each and record in a histogram
	int num_trials = 500;
	BOOST_TEST_MESSAGE("Testing time efficiency of FP & DP on " << num_trials << " instances.");
	double Tmax = 1.0, lb[2*NDOF+1], ub[2*NDOF+1];
	set_bounds(lb,ub,0.01,Tmax);
	joint qact;
	spline_params poly;
	TableTennis tt = TableTennis(true,false);
	arma_rng::set_seed_random();
	//double std_noise = 0.0001;
	//double std_model = 0.3;
	int ball_launch_side;
	int joint_init_pose;
	optim_des racket_params;
	int N = 1000;
	racket_params.Nmax = 1000;
	double time_land_des = 0.8;
	vec2 ball_land_des = {0.0, dist_to_table - 3*table_length/4};
	mat balls_pred = zeros<mat>(6,N);
	vec time_elapsed_init = zeros<vec>(num_trials);
	vec time_elapsed_run_fp = zeros<vec>(num_trials);
	vec time_elapsed_dp_lookup = zeros<vec>(num_trials);
	vec time_elapsed_dp_no_lookup = zeros<vec>(num_trials);
    wall_clock timer;

	for (int n = 0; n < num_trials; n++) { // for each trial
		BOOST_TEST_MESSAGE("Trial: " << n+1);
		ball_launch_side = (randi(1,distr_param(0,2)).at(0));
		joint_init_pose = (randi(1,distr_param(0,2)).at(0));
		init_posture(qact.q,joint_init_pose,false);
		tt.reset_stats();
		tt.set_ball_gun(0.05,ball_launch_side);
		for (int i = 0; i < N; i++) {
			tt.integrate_ball_state(DT);
			balls_pred.col(i) = tt.get_ball_state();
		}
		timer.tic();
		racket_params = calc_racket_strategy(balls_pred,ball_land_des,time_land_des,racket_params);
		FocusedOptim fp = FocusedOptim(qact.q.memptr(),lb,ub);
		fp.set_verbose(false);
		fp.set_des_params(&racket_params);
		fp.update_init_state(qact);
		time_elapsed_init(n) = timer.toc() * 1e3;
		timer.tic();
		fp.run();
		time_elapsed_run_fp(n) = timer.toc() * 1e3;
		DefensiveOptim dp = DefensiveOptim(qact.q.memptr(),lb,ub,true,false);
		dp.set_verbose(false);
		dp.set_des_params(&racket_params);
		dp.update_init_state(qact);
		timer.tic();
		dp.run();
		time_elapsed_dp_no_lookup(n) = timer.toc() * 1e3;
		DefensiveOptim dp2 = DefensiveOptim(qact.q.memptr(),lb,ub,true,true);
		dp2.set_verbose(false);
		dp2.set_des_params(&racket_params);
		dp2.update_init_state(qact);
		timer.tic();
		dp2.run();
		time_elapsed_dp_lookup(n) = timer.toc() * 1e3;
	}
	umat hists = zeros<umat>(15,4);
	hists.col(0) = hist(time_elapsed_init, linspace<vec>(0,9,15));
	hists.col(1) = hist(time_elapsed_run_fp, linspace<vec>(5,75,15));
	hists.col(2) = hist(time_elapsed_dp_no_lookup, linspace<vec>(5,75,15));
	hists.col(3) = hist(time_elapsed_dp_lookup, linspace<vec>(5,75,15));
	hists.save("histograms.txt", raw_ascii);
}
