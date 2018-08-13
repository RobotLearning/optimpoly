'''
Test camera calibration for cameras 0 and 1
Load the pickled dictionary from pixels to 3d ball pos
Perform regression on data and test on still ball data on table
Check out RANSAC as well, is it more robust?
Also estimating projection matrices to compare.
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

pixels = np.zeros((len(ball_locs), 4))
pos3d = np.zeros((len(ball_locs), 3))
for i, tuples in enumerate(ball_locs.values()):
    pixel = tuples[0]
    pos = tuples[1]
    pixels[i, :] = np.array(pixel)
    pos3d[i, :] = np.array(pos)

# pixels = np.arrapos3d(ball_locs.kepos3ds())
pixelsbar = np.hstack((np.ones((pixels.shape[0], 1)), pixels))
# pos3d = np.arrapos3d(ball_locs.values())
sol = linalg.lstsq(pixelsbar, pos3d)
beta = sol[0]
res = sol[1]

# compare with RANSAC
ransac = linear_model.RANSACRegressor()
ransac.fit(pixelsbar, pos3d)
inliers_ran = ransac.inlier_mask_
sol_ran = linalg.lstsq(pixelsbar[inliers_ran, :], pos3d[inliers_ran])
beta_ran = sol_ran[0]
res_ran = sol_ran[1]

# Regressing for projection matrix!!!
X_prj = np.hstack((pos3d, np.ones((pos3d.shape[0], 1))))
y_prj_0 = np.hstack((pixels[:, 0:2], np.ones((pixels.shape[0], 1))))
y_prj_1 = np.hstack((pixels[:, 2:], np.ones((pixels.shape[0], 1))))
sol_prj_mat_0 = linalg.lstsq(X_prj, y_prj_0)
sol_prj_mat_1 = linalg.lstsq(X_prj, y_prj_1)
prj_mat_0 = sol_prj_mat_0[0].T
res_prj_0 = sol_prj_mat_0[1]
prj_mat_1 = sol_prj_mat_1[0].T
res_prj_1 = sol_prj_mat_1[1]
prj_mat = np.vstack((prj_mat_0, prj_mat_1))

# predict 3d ball pos on still balls
# find ball locations for cameras 0 and 1
# check if predictions make sense
# for instance table_length = 2.74 m
img_path = os.environ['HOME'] + '/Dropbox/capture_train/still'
ball_dict = fball.find_balls(
    img_path, ranges=[1, 11], cams=[0, 1], prefix='cam')
pixels_pred = np.array(ball_dict.values())
pixels_pred_bar = np.hstack((np.ones((pixels_pred.shape[0], 1)), pixels_pred))


def test_score(loc_pred):
    ''' output total z-difference squared'''
    o = np.diff(loc_pred[:, -1])
    return np.sum(o*o)


pred_regr_ran_still = np.dot(pixels_pred_bar, beta_ran)
print('train res:', res_ran, 'score:', test_score(pred_regr_ran_still))
pred_regr_still = np.dot(pixels_pred_bar, beta)
print('train_res:', res, 'score:', test_score(pred_regr_still))
# pred_proj_still = linalg.solve(prj_mat[[0, 1, 3, 4], :], pixels_pred.T)[:,0:-1]
import cv2
projPoints0 = pixels_pred[:, 0:2].astype(float).T
projPoints1 = pixels_pred[:, 2:].astype(float).T
points4d = cv2.triangulatePoints(
    prj_mat_0.astype(float), prj_mat_1.astype(float), projPoints0, projPoints1)
pred_proj_still = points4d[0:-1, :].T
for i in range(points4d.shape[1]):
    pred_proj_still[i, :] = pred_proj_still[i, :]/points4d[-1, i]
print('score:', test_score(pred_proj_still))

# compare with old calibration
json_file = os.environ['HOME'] + \
    "/vision/ball_tracking/examples/tracking/server_3d_conf.json"
with open(json_file, 'r') as f:
    old_calib_file = json.load(f)

calibs = old_calib_file["stereo"]["calib"]
calib0 = np.array(calibs[0]['val'])
calib1 = np.array(calibs[1]['val'])
points4d_old = cv2.triangulatePoints(calib0.astype(
    float), calib1.astype(float), projPoints0, projPoints1)
pred_proj_still_old = points4d_old[0:-1, :].T
for i in range(points4d_old.shape[1]):
    pred_proj_still_old[i, :] = pred_proj_still_old[i, :]/points4d_old[-1, i]
print('score:', test_score(pred_proj_still_old))
