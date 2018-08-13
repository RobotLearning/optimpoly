'''
Estimate the P matrix:
Form the A matrix : 2N x 12 matrix
Estimate P with the right singular vector of V in svd(A)
Update with a nonlinear estimation routine
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
import calibrate_nonlin as cal_non


def load_pixels_and_pos():
    ''' Load pixels of both cameras and 3d positions from pickled dictionary'''
    ball_locs = dict()
    img_range = [750, 7380]
    cam_range = [0, 1]
    pickle_file = "python/ball_locs_" + \
                  str(img_range) + "_" + str(cam_range) + ".pickle"
    file_obj = open(pickle_file, 'r')
    ball_locs = pickle.load(file_obj)
    file_obj.close()
    pixels_0 = np.zeros((len(ball_locs), 2))
    pixels_1 = np.zeros((len(ball_locs), 2))
    pos3d = np.zeros((len(ball_locs), 3))
    for i, tuples in enumerate(ball_locs.values()):
        pixel = tuples[0]
        pos = tuples[1]
        pixels_0[i, :] = np.array(pixel[0:2])
        pixels_1[i, :] = np.array(pixel[2:])
        pos3d[i, :] = np.array(pos)
    return pixels_0, pixels_1, pos3d


def estimate_proj_mat_linear(pixels, pos3d):
    ''' Linear estimation of P matrix using SVD decomposition of
    A matrix '''

    N = pixels.shape[0]

    # normalize the images
    mean_pix = np.sum(pixels, axis=0)/N
    d_bar = np.sum(
        np.sqrt((pixels[:, 0]-mean_pix[0])**2 + (pixels[:, 1]-mean_pix[1])**2))
    T = np.zeros((3, 3))
    T[0, 0] = np.sqrt(2)/d_bar
    T[1, 1] = T[0, 0]
    T[2, 2] = 1.0
    T[0, 2] = -np.sqrt(2) * mean_pix[0]/d_bar
    T[1, 2] = -np.sqrt(2) * mean_pix[1]/d_bar

    # normalize the 3d positions
    mean_pos = np.sum(pos3d, axis=0)/N
    D_bar = np.sum(np.sqrt((pos3d[:, 0]-mean_pos[0])**2 +
                           (pos3d[:, 1]-mean_pos[1])**2 + (pos3d[:, 2]-mean_pos[2])**2))
    U = np.zeros((4, 4))
    U[0, 0] = U[1, 1] = U[2, 2] = np.sqrt(3)/D_bar
    U[3, 3] = 1.0
    U[0, 3] = -np.sqrt(3)*mean_pos[0]/D_bar
    U[1, 3] = -np.sqrt(3)*mean_pos[1]/D_bar
    U[2, 3] = -np.sqrt(3)*mean_pos[2]/D_bar

    # form the A matrices
    pixels = np.dot(T, np.vstack((pixels.T, np.ones((1, N)))))
    pos3d = np.dot(U, np.vstack((pos3d.T, np.ones((1, N)))))

    #pixels = np.vstack((pixels.T, np.ones((1, N))))
    #pos3d = np.vstack((pos3d.T, np.ones((1, N))))
    A = np.zeros((2*N, 12))
    for i in range(N):
        a = pos3d[:, i]  # a = np.hstack((pos3d[:, i], 1.0))
        A[2*i, 0:4] = a
        A[2*i, 8:] = -pixels[0, i]*a
        A[2*i+1, 4:8] = a
        A[2*i+1, 8:] = -pixels[1, i]*a

    _, S, Vh = np.linalg.svd(A, full_matrices=True)
    P = Vh[-1, :]
    P = P.reshape((3, 4), order='C')
    # renormalize
    P = np.linalg.solve(T, P.dot(U))
    return P


def test_score(loc_pred):
    ''' output total z-difference squared'''
    o = np.diff(loc_pred[:, -1])
    return np.sum(o*o)


def eval_proj_error(P, pts2d, pts3d):
    ''' Return residual of fitting to pixels given theta parameters
    a.k.a projection error'''

    N = pts3d.shape[0]
    pts4d = np.vstack((pts3d.T, np.ones((1, N))))
    proj_pts = np.dot(P, pts4d)
    difs = pts2d.T - proj_pts[0:-1, :]
    res = np.sum(difs*difs)
    print('residual before nonlin opt:', res)


def eval_on_still_balls(P0, P1):
    ''' Evaluate camera models by triangulating to predict still balls'''

    # predict 3d ball pos on still balls
    # find ball locations for cameras 0 and 1
    # check if predictions make sense
    # for instance table_length = 2.74 m
    img_path = os.environ['HOME'] + '/Dropbox/capture_train/still'
    ball_dict = fball.find_balls(
        img_path, ranges=[1, 11], cams=[0, 1], prefix='cam')
    pixels = np.array(ball_dict.values())

    points4d = cv2.triangulatePoints(P0.astype(
        float), P1.astype(float), pixels[:, 0:2].astype(float).T, pixels[:, 2:].astype(float).T)
    pred_proj_still = points4d[0:-1, :].T
    pred_proj_still_norm = np.zeros(pred_proj_still.shape)
    for i in range(points4d.shape[1]):
        pred_proj_still_norm[i, :] = pred_proj_still[i, :]/points4d[-1, i]
    print('pred 3d points:')
    print(pred_proj_still_norm)
    print('score:', test_score(pred_proj_still_norm))


pixels_0, pixels_1, pos3d = load_pixels_and_pos()
P0 = estimate_proj_mat_linear(pixels_0, pos3d)
P1 = estimate_proj_mat_linear(pixels_1, pos3d)
eval_on_still_balls(P0, P1)
eval_proj_error(P0, pixels_0, pos3d)
eval_proj_error(P1, pixels_1, pos3d)

''' PREPARE FOR A NONLINEAR LEAST SQUARES BASED CALIB UPDATE'''
out0 = cv2.decomposeProjectionMatrix(P0)
out1 = cv2.decomposeProjectionMatrix(P1)
cam_mat_0 = out0[0]
trans_vec_0 = out0[2].reshape((4,))
trans_vec_0 = trans_vec_0[0:-1]/trans_vec_0[-1]
cam_mat_1 = out1[0]
trans_vec_1 = out1[2].reshape((4,))
trans_vec_1 = trans_vec_1[0:-1]/trans_vec_1[-1]
intrinsic_dict_0 = {'fx': cam_mat_0[0, 0], 'fy': cam_mat_0[1, 1],
                    'shear': cam_mat_0[0, 1], 'u0': cam_mat_0[0, 2], 'v0': cam_mat_0[1, 2]}
intrinsic_dict_1 = {'fx': cam_mat_1[0, 0], 'fy': cam_mat_1[1, 1],
                    'shear': cam_mat_1[0, 1], 'u0': cam_mat_1[0, 2], 'v0': cam_mat_1[1, 2]}
extrinsic_dict_0 = {'euler_angles': out0[-1]
                    * 2*np.pi/360, 'trans_vector': trans_vec_0}
extrinsic_dict_1 = {'euler_angles': out1[-1]
                    * 2*np.pi/360, 'trans_vector': trans_vec_1}

params_0 = cal_non.est_calib_params_nonlin(distortion_dict=None,
                                           intrinsic_dict=intrinsic_dict_0,
                                           extrinsic_dict=extrinsic_dict_0,
                                           pts3d=pos3d.T,
                                           pts2d=pixels_0.T,
                                           num_iter_max=1000,
                                           debug=False)
params_1 = cal_non.est_calib_params_nonlin(distortion_dict=None,
                                           intrinsic_dict=intrinsic_dict_1,
                                           extrinsic_dict=extrinsic_dict_1,
                                           pts3d=pos3d.T,
                                           pts2d=pixels_1.T,
                                           num_iter_max=1000,
                                           debug=False)

# form P0_new and P1_new
# undistort the pixels
'''
P0_new = params_0['proj']
P1_new = params_1['proj']


dist_coeffs_0 = np.array(params_0['dist'].values()).astype(np.float32)
dist_coeffs_1 = np.array(params_1['dist'].values()).astype(np.float32)
pixels_0_undistort = np.zeros(pixels_0.shape)
pixels_1_undistort = np.zeros(pixels_0.shape)


# dist_coeffs_0 = np.ones((1, 8), dtype=np.float32)
# dist_coeffs_1 = np.ones((1, 8), dtype=np.float32)
pixels_0_undistort = cv2.undistortPoints(pixels_0[:, np.newaxis, :].astype(
    np.float32), cam_mat_0.astype(np.float32), dist_coeffs_0)
pixels_1_undistort = cv2.undistortPoints(pixels_1[:, np.newaxis, :].astype(
    np.float32), cam_mat_1.astype(np.float32), dist_coeffs_1)

#eval_proj_error(P0_new, pixels_0, pos3d)
#eval_proj_error(P1_new, pixels_1, pos3d)
eval_on_still_balls(P0_new, P1_new)
'''
