'''
load data
cut it and process
plot an example
'''

import numpy as np
import argparse
import barrett_wam_kinematics as wam
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from mpl_toolkits.mplot3d import Axes3D

parser = argparse.ArgumentParser(
    description='Load saved joint and ball data from demonstrations. Process it and plot.')
parser.add_argument('--joint_file', help='joint file')
parser.add_argument('--ball_file', help='ball file')
parser.add_argument(
    '--num_examples', help='number of demonstrations', type=int)
parser.add_argument('--plot_example', help='plot a specific example', type=int)
parser.add_argument(
    '--plot_3d', help='also plot the cartesian values', action="store_true")
parser.add_argument(
    '--add_time', help='joints and ball files includes absolutime time')
parser.add_argument(
    '--smooth', help='smoothing factor of splines while plotting')
args = parser.parse_args()
assert (args.plot_example < 
        args.num_examples), "example to plot must be less than num of examples"

joint_file = args.joint_file  # './data/10.6.18/joints.txt'
ball_file = args.ball_file
num_examples = args.num_examples  # 21  # first dataset = 20, second = 21
M = np.fromfile(joint_file, sep=" ")
ndof = 7
dt = 0.002
if args.add_time == 0:
    N = M.size/ndof
    q = np.reshape(M, [N, ndof])
    t_joints = dt * np.linspace(1, N, N)
else:
    N = M.size/(ndof+1)
    M = np.reshape(M, [N, ndof+1])
    q = M[:,1:]
    t_joints = M[:,0]
    t_joints = 0.001 * t_joints #in miliseconds

t_min = t_joints[0]
if args.ball_file:
    B = np.fromfile(ball_file, sep=" ")
    if args.add_time == 0:
        N_balls = B.size/3
        balls = np.reshape(B, [N_balls, 3])
        t_balls = dt * np.linspace(1, N, N)
    else:
        N_balls = B.size/4
        B = np.reshape(B, [N_balls, 4])
        balls = B[:,1:]
        t_balls = B[:,0]
        t_balls = 0.001 * t_balls
        t_min = min(t_min, t_balls[0])
        t_balls = t_balls - t_min
        t_joints = t_joints - t_min

# quick plotting for ball
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(balls[:,0], balls[:,1], balls[:,2], c="b")
plt.show()


def detect_movements(x, num_examples, dt):
    '''
    Detect movements by checking for maximum velocity
    x can be ball data or joint data
    '''

    # Find the prominent fast moving segments
    xd = np.diff(x, 1, 0)/dt
    vels = np.sqrt(np.sum(xd*xd, -1))
    sorted_vels = np.sort(vels)
    sorted_vels = sorted_vels[::-1]
    idx_vels = np.argsort(vels)
    idx_vels = idx_vels[::-1]

    idx_clusters = np.array([idx_vels[0]])
    idx = 1
    thresh = 500  # 1 sec diff min betw demonstr.
    while idx_clusters.size < num_examples:
        diff_max = min(abs(idx_vels[idx] - idx_clusters))
        if diff_max > thresh:
            idx_clusters = np.insert(idx_clusters, 0, idx_vels[idx]+1)
        idx = idx+1

    # sort the points and find points of low velocity around them
    clusters = np.sort(idx_clusters)
    low_vel_idxs = np.zeros((2, num_examples))
    low_vel_thresh = 1e-3
    max_duration_movement = 1.0  # seconds
    idx_max_move = max_duration_movement/dt
    #print clusters
    for i in range(num_examples):
        # find the first index below cluster high vel idx where
        # the vel drops below thresh
        idx_low_pt = clusters[i]
        iter = 0
        while vels[idx_low_pt] > low_vel_thresh and iter < idx_max_move/2:
            idx_low_pt = idx_low_pt - 1
            iter = iter + 1
            low_vel_idxs[0, i] = idx_low_pt

        # find the first index above cluster idx
        idx_high_pt = clusters[i]
        iter = 0
        while vels[idx_high_pt] > low_vel_thresh and iter < idx_max_move/2:
            idx_high_pt = idx_high_pt + 1
            iter = iter + 1
            low_vel_idxs[1, i] = idx_high_pt

    return low_vel_idxs


def plot_examples(examples, t, x, idx_movements, smooth_fact=0.01, dim=7, plot_3d=False):

    num_examples = examples.size
    for i in range(num_examples):
        f, axs = plt.subplots(7, 1, sharex=False)
        # print examples[i]
        idx_plot = np.arange(
            start=idx_movements[0, examples[i]],
            stop=idx_movements[1, examples[i]]+1, step=1, dtype=np.int32)
        print idx_plot
        x_plot = x[idx_plot, :]
        t_plot = t[idx_plot]
        for j in range(dim):
            spl = UnivariateSpline(
                t_plot, x_plot[:, j], w=None, k=3, s=smooth_fact)
            x_smooth = spl(t_plot)
            axs[j].plot(t_plot, x_plot[:, j])
            axs[j].plot(t_plot, x_smooth)
        # KINEMATICS PLOT
        if plot_3d:
            x_plot = np.zeros((3, idx_plot.size))
            for idx, val in enumerate(idx_plot):
                As = wam.barrett_wam_kinematics(np.transpose(x_plot[idx, :]))
                x_plot[:, idx] = np.transpose(As[-1])[-1][0:3]
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(x_plot[0, :], x_plot[1, :], x_plot[2, :], c="b")
        plt.show()

'''
# plot one of the motions
#
#GOOD MOTIONS (21 TOTAL)
#[2,4,6,9,10,12,15,16,18,19,20,21]-1 (zero indexing)

idx_joint_movements = detect_movements(q, num_examples, dt)
if args.ball_file:
    idx_ball_movements = detect_movements(balls, num_examples, dt=1/180.0)
#print idx_movements

if args.smooth:
    smooth = args.smooth
else:
    smooth = 0.01

#examples = np.array([8])
if args.plot_example:
    examples = np.array([args.plot_example])
    if args.plot_3d:
        plot_examples(examples, t_joints, q, idx_joint_movements, smooth_fact=smooth, 
                      plot_3d=True, dim=7)
        if args.ball_file:
            plot_examples(examples, t_balls, balls, dim=3,
                          idx_movements=idx_ball_movements, smooth_fact=smooth)
    else:
        plot_examples(examples, t_joints, q, idx_joint_movements, smooth_fact=smooth)
'''
