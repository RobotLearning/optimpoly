'''
Estimate the P matrix:
Form the A matrix : 2N x 12 matrix
Estimate P with the right singular vector of V in svd(A)
Update with a nonlinear estimation routine
Predict on still balls on the table
'''

import pickle
import numpy as np
import scipy.linalg as linalg
import json
import os
import sys
#from sklearn import linear_model
import find_balls as fball
import cv2
import calibrate_nonlin as cal_non
sys.path.append('./python')


def load_pixels_and_pos(dataset, img_range, cam_range):
    ''' Load pixels of both cameras and 3d positions from pickled dictionary'''
    ball_locs = dict()
    pickle_file = "python/ball_locs_" + dataset + '_' + \
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


def normalize_pixels_and_pos3d(pixels, pos3d):
    ''' normalize the matrices before SVD '''

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
    return pixels, pos3d, T, U


def estimate_proj_mat_linear(pixels, pos3d):
    ''' Linear estimation of P matrix using SVD decomposition of
    A matrix '''

    N = pixels.shape[0]
    pixels, pos3d, T, U = normalize_pixels_and_pos3d(pixels, pos3d)
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


def test_score_over_table(loc_pred):
    '''
    test prediction accuracy over table by calculating
    z- score: total deviation
    y score: average distance of points on robot court [0-8 inclusive]
             should be roughly table_length/4
    x score: balls are on the edge always so, [-x,0,x] mesh

    Geometry of placed balls
    [0 1 2
     5 4 3
     6 7 8 (close to net)
     ----- [net]
      10 9 (close to net)
    '''
    table_length = 2.76
    table_width = 1.525
    x_edge = table_width/2.0
    x_center = 0.0
    center_balls = [1, 4, 7]
    left_balls = [0, 5, 6]
    right_balls = [2, 3, 8]

    xdifs = loc_pred[center_balls, 0] - x_center + \
        loc_pred[left_balls, 0] - (-x_edge) + \
        loc_pred[right_balls, 0] - x_edge
    ydifs = np.sum(
        np.abs(loc_pred[[0, 1, 2], 1] - loc_pred[[5, 4, 3], 1]) - table_length/4.0 +
        np.abs(loc_pred[[5, 4, 3], 1] - loc_pred[[6, 7, 8], 1]) - table_length/4.0)
    zdifs = np.diff(loc_pred[:, -1])
    return np.sum(zdifs*zdifs) + np.sum(xdifs*xdifs) + np.sum(ydifs*ydifs)


def eval_proj_error(P, pts2d, pts3d):
    ''' Return residual of fitting to pixels given theta parameters
    a.k.a projection error'''

    N = pts3d.shape[0]
    pts4d = np.vstack((pts3d.T, np.ones((1, N))))
    proj_pts = np.dot(P, pts4d)
    difs = pts2d.T - proj_pts[0:-1, :]
    res = np.sum(difs*difs)
    print('residual:', res)


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
    print('pixels for predicting 3d points:')
    print(pixels)
    P0 = P0.astype(float)
    P1 = P1.astype(float)
    pixels = pixels.astype(float)
    points4d = cv2.triangulatePoints(
        P0, P1, pixels[:, 0:2].T, pixels[:, 2:].T)
    # normalize
    points3d = points4d[0:-1, :] / points4d[-1, :]

    # LINEAR ALGEBRAIC TRIANGULATION
    # is exactly the same as cv2.triangulatePoints
    '''
    A = np.zeros((4, 4))
    points4d_mine = np.zeros(points4d.shape)
    for i in range(11):
        A[0, :] = pixels[i, 0]*P0[-1, :] - P0[0, :]
        A[1, :] = pixels[i, 1]*P0[-1, :] - P0[1, :]
        A[2, :] = pixels[i, 2]*P1[-1, :] - P1[0, :]
        A[3, :] = pixels[i, 3]*P1[-1, :] - P1[1, :]
        U, S, Vh = np.linalg.svd(A, full_matrices=True)
        points4d_mine[:, i] = Vh[-1, :]
    points3d_mine = points4d_mine[0:-1, :] / points4d_mine[-1, :]
    print(points3d_mine.T)
    '''

    print('pred 3d points:')
    print(points3d.T)
    print('score over table:', test_score_over_table(points3d.T))


def decompose_proj_mat(P):
    ''' Decompose projection matrix into calib mat A, rotation R and translation t'''
    '''
    # BUGGY METHOD
    B = P[:, 0:-1]
    b = P[:, -1]
    K = B.dot(B.T)
    scale = K[2, 2]
    K = K / K[2, 2]  # normalize
    A = np.zeros((3, 3))
    A[0, 2] = K[0, 2]  # u0
    A[1, 2] = K[1, 2]  # v0
    A[1, 1] = np.sqrt(K[1, 1] - A[1, 2]**2)  # fy
    A[0, 1] = (K[1, 0] - (A[0, 2]*A[1, 2]))/A[1, 1]  # shear
    A[0, 0] = np.sqrt(K[0, 0] - A[0, 2]**2 - A[0, 1]**2)  # fx
    A[2, 2] = 1.0
    R = np.linalg.solve(A, B)
    t = np.linalg.solve(A, b)
    return A, R, t, scale
    '''

    # use instead RQ transform
    K, R = linalg.rq(P[:, :-1])
    t = linalg.solve(K, P[:, -1])
    # scale = K[-1, 1]
    return K, R, t


def form_rot_matrix(euler_angles):
    ''' Form the rotation matrix from 3 euler angles
    for some reson we have to also take transpose!!'''
    mat0 = np.eye(3)
    mat0[1, 1] = np.cos(euler_angles[0])
    mat0[1, 2] = np.sin(euler_angles[0])
    mat0[2, 1] = -mat0[1, 2]
    mat0[2, 2] = mat0[1, 1]

    mat1 = np.eye(3)
    mat1[0, 0] = np.cos(euler_angles[1])
    mat1[0, 2] = -np.sin(euler_angles[1])
    mat1[2, 0] = -mat1[0, 2]
    mat1[2, 2] = mat1[0, 0]

    mat2 = np.eye(3)
    mat2[0, 0] = np.cos(euler_angles[2])
    mat2[0, 1] = np.sin(euler_angles[2])
    mat2[1, 0] = -mat2[0, 1]
    mat2[1, 1] = mat2[0, 0]

    return np.dot(mat0, np.dot(mat1, mat2)).T


def form_extrinsic_mat(extrinsic_dict):
    ''' form the matrix composed of rotation and translation (scaling is 1)'''
    euler_angles = extrinsic_dict['euler_angles']
    trans = extrinsic_dict['trans_vector'][np.newaxis].T
    rot_mat = form_rot_matrix(euler_angles)
    return np.hstack((rot_mat, trans))


def form_intrinsic_mat(intrinsic_dict):
    ''' forming the camera (intrinsic) matrix, and scaling'''
    mat = np.zeros((3, 3))
    mat[0, 0] = intrinsic_dict['fx']
    mat[0, 1] = intrinsic_dict['shear']
    mat[0, 2] = intrinsic_dict['u0']
    mat[1, 1] = intrinsic_dict['fy']
    mat[1, 2] = intrinsic_dict['v0']
    mat[2, 2] = intrinsic_dict['scale']
    return mat


pixels_0_red, pixels_1_red, pos3d_red = load_pixels_and_pos(
    dataset='red', img_range=(1510, 6080), cam_range=(0, 1))
pixels_0_black, pixels_1_black, pos3d_black = load_pixels_and_pos(
    dataset='black', img_range=(750, 7380), cam_range=(0, 1))

pixels_0 = pixels_0_red
pixels_1 = pixels_1_red
pos3d = pos3d_red
P0 = estimate_proj_mat_linear(pixels_0, pos3d)
P1 = estimate_proj_mat_linear(pixels_1, pos3d)
eval_on_still_balls(P0, P1)
eval_proj_error(P0, pixels_0_red, pos3d_red)
eval_proj_error(P1, pixels_1_red, pos3d_red)
eval_proj_error(P0, pixels_0_black, pos3d_black)
eval_proj_error(P1, pixels_1_black, pos3d_black)
# evaluate on the other (black/red) dataset


# PREPARE FOR A NONLINEAR LEAST SQUARES BASED CALIB UPDATE

A0, R0, t0 = decompose_proj_mat(P0)
A1, R1, t1 = decompose_proj_mat(P1)

'''
# USING SOLVEPNP TO REFINE DOESNT HELP
R0_vec = cv2.Rodrigues(R0)[0]
t0_new = t0.copy()
_ = cv2.solvePnPRansac(pos3d, pixels_0, cameraMatrix=A0, rvec=R0_vec, tvec=t0_new,
                       distCoeffs=np.zeros((4,
                                            1)),
                       flags=cv2.SOLVEPNP_EPNP,  # cv2.SOLVEPNP_ITERATIVE,
                       useExtrinsicGuess=True)
R0_new, _ = cv2.Rodrigues(R0_vec)
P0_new = A0.dot(np.hstack((R0_new, t0_new.reshape((3, 1)))))
eval_proj_error(P0_new, pixels_0, pos3d)

R1_vec = cv2.Rodrigues(R1)[0]
t1_new = t1.copy()
_ = cv2.solvePnPRansac(pos3d, pixels_1, cameraMatrix=A1, rvec=R1_vec, tvec=t1_new,
                       distCoeffs=np.zeros((4,
                                            1)),
                       flags=cv2.SOLVEPNP_EPNP,  # cv2.SOLVEPNP_ITERATIVE,
                       useExtrinsicGuess=True)
R1_new, _ = cv2.Rodrigues(R1_vec)
P1_new = A1.dot(np.hstack((R1_new, t1_new.reshape((3, 1)))))
eval_proj_error(P1_new, pixels_1, pos3d)

eval_on_still_balls(P0_new, P1_new)
'''

# CUSTOM NONLINEAR LS CODE
'''
cam_mat_0 = A0.copy()
cam_mat_1 = A1.copy()
out0 = cv2.decomposeProjectionMatrix(P0)
out1 = cv2.decomposeProjectionMatrix(P1)
intrinsic_dict_0 = {'fx': cam_mat_0[0, 0], 'fy': cam_mat_0[1, 1],
                    'shear': cam_mat_0[0, 1], 'u0': cam_mat_0[0, 2],
                    'v0': cam_mat_0[1, 2], 'scale': cam_mat_0[2, 2]}
intrinsic_dict_1 = {'fx': cam_mat_1[0, 0], 'fy': cam_mat_1[1, 1], 'shear':
                    cam_mat_1[0, 1], 'u0': cam_mat_1[0, 2], 'v0':
                    cam_mat_1[1, 2], 'scale': cam_mat_1[2, 2]}
extrinsic_dict_0 = {'euler_angles': out0[-1] * 2*np.pi/360,
                    'trans_vector': t0}
extrinsic_dict_1 = {'euler_angles': out1[-1] * 2*np.pi/360,
                    'trans_vector': t1}

params_0 = cal_non.est_calib_params_nonlin(distortion_dict=None,
                                           intrinsic_dict=intrinsic_dict_0,
                                           extrinsic_dict=extrinsic_dict_0,
                                           pts3d=pos3d.T,
                                           pts2d=pixels_0.T,
                                           num_iter_max=10000,
                                           debug=False)
params_1 = cal_non.est_calib_params_nonlin(distortion_dict=None,
                                           intrinsic_dict=intrinsic_dict_1,
                                           extrinsic_dict=extrinsic_dict_1,
                                           pts3d=pos3d.T,
                                           pts2d=pixels_1.T,
                                           num_iter_max=10000,
                                           debug=False)

# form P0_new and P1_new

A0_new = form_intrinsic_mat(params_0['intrinsic'])
E0_new = form_extrinsic_mat(params_0['extrinsic'])
P0_new = A0_new.dot(E0_new)

A1_new = form_intrinsic_mat(params_1['intrinsic'])
E1_new = form_extrinsic_mat(params_1['extrinsic'])
P1_new = A1_new.dot(E1_new)

eval_proj_error(P0_new, pixels_0, pos3d)
eval_proj_error(P1_new, pixels_1, pos3d)
eval_on_still_balls(P0_new, P1_new)
'''

# undistort the pixels

'''
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

'''

# compare with old calibration
'''
print('Comparing with old calibration...')
json_file = os.environ['HOME'] + \
    "/table-tennis/json/server_3d_conf_ping.json"
with open(json_file, 'r') as f:
    old_calib_file = json.load(f)

calibs = old_calib_file["stereo"]["calib"]
calib0 = np.array(calibs[0]['val'])
calib1 = np.array(calibs[1]['val'])

# the robot loc changed since we calibrated last time
dx = np.array([0.0, 0.0, 0.0])
R0 = calib0[:, 0:-1]
R1 = calib1[:, 0:-1]
dt0 = -np.dot(R0, dx)
dt1 = -np.dot(R1, dx)
calib0[:, -1] = calib0[:, -1] + dt0
calib1[:, -1] = calib1[:, -1] + dt1

eval_proj_error(calib0, pixels_0, pos3d)
eval_proj_error(calib1, pixels_1, pos3d)
eval_on_still_balls(calib0, calib1)
'''
