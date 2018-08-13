'''
DOES NOT WORK!! REQUIRES CHESSBOARD PATTERN!!
Test nonlinear least squares based calibration
with opencv calibrateCamera function
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
import cv2

ball_locs = dict()
img_range = [750, 7380]
cam_range = [0, 1]
pickle_file = "python/ball_locs_" + \
    str(img_range) + "_" + str(cam_range) + ".pickle"
file_obj = open(pickle_file, 'r')
ball_locs = pickle.load(file_obj)
file_obj.close()

X = np.zeros((len(ball_locs), 4))
y = np.zeros((len(ball_locs), 3))
for i, tuples in enumerate(ball_locs.values()):
    pixels = tuples[0]
    pos = tuples[1]
    X[i, :] = np.array(pixels)
    y[i, :] = np.array(pos)

objpoints = []
imgpoints_0 = []
imgpoints_1 = []

for i in range(len(ball_locs)):
    objpoints.append(y[i, :][np.newaxis].astype('float32'))
    imgpoints_0.append(X[i, 0:2][np.newaxis].astype('float32'))
    imgpoints_1.append(X[i, 2:][np.newaxis].astype('float32'))

shape = (659, 494)
flags = (cv2.CALIB_SAME_FOCAL_LENGTH + cv2.CALIB_ZERO_TANGENT_DIST)
rot, trans, ess, fundam = cv2.stereoCalibrate(
    [np.array(objpoints)], [np.array(imgpoints_0)], [np.array(imgpoints_1)], imageSize=shape, cameraMatrix1=None, distCoeffs1=None, distCoeffs2=None, cameraMatrix2=None, flags=flags)
