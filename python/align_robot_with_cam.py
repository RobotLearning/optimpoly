'''
1. Detect start and finish of movement in joint space
2. check alignment with camera data
camera 0: robot starts moving around 750 and stops around 7380
for camera 1 it's the same
3. for every camera frame in between run kinematics to get 3d position of racket
and also orientation
ball is 14 cm behind racket and 5 cm on top
for black racket side on top 5 cm should be subtracted
'''
import sys
import os
sys.path.append("python/")
sys.path.append(os.environ['HOME'] + "/vision/Vision_Experiments/sebastian/")
import numpy as np
from scipy import signal
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import barrett_wam_kinematics as wam
import racket_calc as racket
import tensorflow as tf
import ball_utils
import pickle
import find_balls as fball


joint_file = './data/29.7.18/joints_black.txt'
M = np.fromfile(joint_file, sep=" ")
ndof = 7
N = M.size/ndof
q = np.reshape(M, [N, ndof])
dt = 0.002
t = dt * np.linspace(1, N, N)

# keep only signal from 20 to 60 seconds
idx_cut = np.logical_and(t >= 20, t <= 60)
tcut = t[idx_cut]
q = q[idx_cut]
t = tcut
# t = t[t>=10]
# q = q[t>=10]

# smoothen the signal
b, a = signal.butter(2, 0.1)
qfilt = np.zeros(q.shape)
for i in range(ndof):
    qfilt[:, i] = signal.filtfilt(b, a, q[:, i], padlen=50)


# plot the first joint with smoothened values
f, axs = plt.subplots(ndof, 1, sharex=False)
for i in range(ndof):
    axs[i].plot(t, q[:, i], 'b-', label='input')
    axs[i].plot(t, qfilt[:, i], 'r-', label='filt')
    plt.legend(loc='best')
# plt.show()

# get velocity
qd = np.diff(qfilt, 1, 0)/dt
vels = np.sqrt(np.sum(qd*qd, -1))
# find first occurrence of velocity exceeding threshold
thresh = 0.1
idx_start = np.nonzero(vels > thresh)[0][0]
idx_end = np.nonzero(vels > thresh)[0][-1]
idx_slack = 1
t_move = t[idx_start-idx_slack:idx_end+idx_slack]
q_move = qfilt[idx_start-idx_slack:idx_end+idx_slack, :]

# plot the first joint with smoothened values
'''
f2, axs2 = plt.subplots(ndof, 1, sharex=False)
for i in range(ndof):
    axs2[i].plot(t_move, q_move[:, i], 'b-', label='input')
    # axs[i].plot(t,qfilt[:,i],'r-',label='filt')
    plt.legend(loc='best')
plt.show()
'''

# get ball positions
# subtract 14 cm along racket slide and add 5 cm along racket normal
x_ball = np.zeros((3, t_move.size))
for idx, tval in enumerate(t_move):
    As = wam.barrett_wam_kinematics(np.transpose(q_move[idx, :]))
    R = As[-1, 0:3, 0:3]
    x_racket = As[-1, 0:3, -1]
    quat = racket.rot2Quat(R)
    orient = racket.calcRacketOrientation(quat)
    R = np.squeeze(racket.quat2Rot(orient))
    x_ball[:, idx] = x_racket + 0.14 * R[:, 1] - 0.05 * R[:, 2]

# load frames and detect centre of balls
frame_path = os.environ['HOME'] + \
    "/Desktop/capture_okan_2907118/movies/robot/ping/black"
img_range = [750, 7380]
cam_range = [0, 1]
ball_dict = fball.find_balls(img_path=frame_path,
                             ranges=img_range,
                             cams=cam_range, debug=False)
img_list = ball_dict.keys()
ball_locs = dict()
# interpolate between 7380 and 750th frames to get 3d positions for that particular frame
for img_num, pixels_ball in ball_dict.iteritems():
    i = img_num - img_list[0]
    t_img = t_move[0] + (i/len(ball_dict)) * (t_move[-1] - t_move[0])
    x_ball_img = []
    for j in range(3):
        x_ball_img.append(np.interp(t_img, xp=t_move, fp=x_ball[j, :]))
    ball_locs[tuple(pixels_ball)] = x_ball_img

pickle_file = "ball_locs_" + str(img_range) + "_" + str(cam_range) + ".pickle"
file_obj = open(pickle_file, 'wb')
pickle.dump(ball_locs, file_obj)
file_obj.close()
