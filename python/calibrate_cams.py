'''
Calibrate cameras 0 and 1
Load the pickled dictionary from pixels to 3d ball pos
Perform regression on half data and validate on other half
Check out RANSAC as well, is it more robust?
'''
import pickle
import numpy as np
import scipy.linalg as linalg
import json
import os
import sys
sys.path.append('./python')
from sklearn import linear_model
import find_balls as fball

ball_locs = dict()
img_range = [750, 7380]
cam_range = [0, 1]
pickle_file = "python/ball_locs_" + \
    str(img_range) + "_" + str(cam_range) + ".pickle"
file_obj = open(pickle_file, 'r')
ball_locs = pickle.load(file_obj)
file_obj.close()

'''
ind = 0
for key, val in ball_locs.iteritems():
    if ind < 10:
        print(key, val)
    else:
        break
    ind = ind + 1
'''
X = np.zeros((len(ball_locs), 4))
y = np.zeros((len(ball_locs), 3))
for i, tuples in enumerate(ball_locs.values()):
    pixels = tuples[0]
    pos = tuples[1]
    X[i, :] = np.array(pixels)
    y[i, :] = np.array(pos)

#X = np.array(ball_locs.keys())
Xbar = np.hstack((np.ones((X.shape[0], 1)), X))
#y = np.array(ball_locs.values())
sol = linalg.lstsq(Xbar, y)
beta = sol[0]
res = sol[1]
# regularized soln
# sol_reg = np.linalg.lstsq(Xbar.T.dot(Xbar) + 4e-2 *
#                          np.identity(Xbar.shape[1]), Xbar.T.dot(y))
#beta_reg = sol_reg[0]

# compare with RANSAC
ransac = linear_model.RANSACRegressor()
ransac.fit(Xbar, y)
inliers_ran = ransac.inlier_mask_
sol_ran = linalg.lstsq(Xbar[inliers_ran, :], y[inliers_ran])
beta_ran = sol_ran[0]
res_ran = sol_ran[1]

# predict 3d ball pos on still balls
# find ball locations for cameras 0 and 1
# check if predictions make sense
# for instance table_length = 2.74 m
img_path = os.environ['HOME'] + '/Dropbox/capture_train/still'
ball_dict = fball.find_balls(
    img_path, ranges=[1, 11], cams=[0, 1], prefix='cam')
X_pred = np.array(ball_dict.values())
X_pred_bar = np.hstack((np.ones((X_pred.shape[0], 1)), X_pred))

loc_pred_ran_still = np.dot(X_pred_bar, beta_ran)
loc_pred_still = np.dot(X_pred_bar, beta)
#loc_pred_still_reg = np.dot(X_pred_bar, beta_reg)
# compare with old calibration
'''
json_file = os.environ['HOME'] + \
    "/vision/ball_tracking/server_3d_conf_ping.json"
with open(json_file, 'r') as f:
    old_calib_file = json.load(f)

calibs = old_calib_file["stereo"]["calib"]
calib0 = np.array(calibs[0]['val'])
calib1 = np.array(calibs[1]['val'])
res_old_0 = calib0.dot(X.transpose()) - y.transpose()
res_old_0 = np.sum(res_old_0 * res_old_0, axis=1)
res_old_1 = calib1.dot(X.transpose()) - y.transpose()
res_old_1 = np.sum(res_old_1 * res_old_1, axis=1)
'''
