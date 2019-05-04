static bool fuse_blobs(const vec3 &blob1, const vec3 &blob3,
                       const bool &status1, const bool &status3,
                       const vec6 &ball_est, vec3 &obs);
static bool check_blob_validity(const vec3 &blob, const bool &status);

/*
 * Test outlier prediction with Extended Kalman Filter
 *
 * Here we want to test/debug the Kalman Filter
 * that is to be used realtime with our vision system
 * for ball tracking
 *
 * Things to watch out for:
 * 1. It should reset correctly everytime ball is launched from ball gun
 * 2. Outliers should be cleaned properly
 * 3. New balls should be updated
 *
 */
void test_outlier_detection() {

  BOOST_TEST_MESSAGE("Testing filtering on REAL BALL DATA!");
  vec3 blob1, blob3, obs;
  bool status1, status3;
  mat real_ball_data;
  bool outlier_detection, predict_with_spin;
  std::string home = std::getenv("HOME");
  std::string trial;
  try {
    cout << "Which trial data to filter?\n";
    cin >> trial;
    cout << "Outlier detection? " << endl;
    cin >> outlier_detection;
    cout << "Predict with spin?" << endl;
    cin >> predict_with_spin;
    cout << "Filtering trial: " << trial << endl;
    real_ball_data.load(home + "/Dropbox/data/real_ball_data_0317/balls_" +
                            trial + ".txt",
                        raw_ascii);
  } catch (const char *exception) {
    cout << "Problem accessing/finding real ball data on Dropbox!" << endl;
  }
  mat ball_states = zeros<mat>(real_ball_data.n_rows, 6);
  player_flags flags;
  flags.outlier_detection = outlier_detection;
  flags.spin = predict_with_spin;
  flags.verbosity = 3;
  flags.min_obs = 12;
  flags.var_model = 0.1;
  flags.var_noise = 0.001;
  flags.detach = true;
  Player cp = Player(zeros<vec>(7), flags);
  for (unsigned i = 0; i < real_ball_data.n_rows; i++) {
    status1 = real_ball_data(i, 1);
    blob1 = real_ball_data(i, span(2, 4)).t();
    status3 = real_ball_data(i, 6);
    blob3 = real_ball_data(i, span(7, 9)).t();
    fuse_blobs(blob1, blob3, status1, status3, ball_states.row(i).t(),
               obs); // get observation
    ball_states.row(i) = cp.filt_ball_state(obs).t();
    usleep(2000);
  }
  ball_states.save(home + "/Dropbox/data/real_ball_data_0317/balls_" + trial +
                       "_filtered.txt",
                   raw_ascii);
}

/*
 *
 * Fusing the blobs
 * If both blobs are valid blob3 is preferred
 * Only updates if the blobs are valid, i.e. not obvious outliers
 *
 */
static bool fuse_blobs(const vec3 &blob1, const vec3 &blob3,
                       const bool &status1, const bool &status3,
                       const vec6 &ball_est, vec3 &obs) {

  bool status = false;
  bool filter_init = ball_est(DY) > 1.0;

  // if ball is detected reliably
  // Here we hope to avoid outliers and prefer the blob3 over blob1
  if (check_blob_validity(blob3, status3) ||
      (check_blob_validity(blob1, status1) && filter_init)) {
    status = true;
    if (status3 && status1 && filter_init)
      obs = (blob3 + blob1) / 2;
    else
      obs = (status3) ? blob3 : blob1;
  }
  return status;
}

/*
 *
 * Checks for the validity of blob ball data using obvious table tennis checks.
 * Returns TRUE if valid.
 *
 * Does not use uncertainty estimates to assess validity
 * so do not rely on it as the sole source of outlier detection!
 *
 */
static bool check_blob_validity(const vec3 &blob, const bool &status) {

  bool valid;
  static vec3 last_blob = blob;
  static double zMax = 0.5;
  static double zMin = floor_level - table_height;
  static double xMax = table_width / 2.0;
  static double yMax = 0.5;
  static double yMin = dist_to_table - table_length - 0.5;
  static double yCenter = dist_to_table - table_length / 2.0;

  if (!status) {
    valid = false;
  } else if (blob(Z) > zMax) {
    valid = false;
  } else if (blob(Y) > yMax || blob(Y) < yMin) {
    valid = false;
  } else if (last_blob(Y) < yCenter && blob(Y) > dist_to_table) {
    valid = false;
  }
  // on the table blob should not appear under the table
  else if (fabs(blob(X)) < xMax &&
           fabs(blob(Y) - yCenter) < table_length / 2.0 && blob(Z) < zMin) {
    valid = false;
  } else {
    valid = true;
  }
  last_blob = blob;
  return valid;
}
