'''
load data
cut it and process
plot an example
'''
import sys
import os
import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from mpl_toolkits.mplot3d import Axes3D
sys.path.append("python/")
import barrett_wam_kinematics as wam
import racket_calc as racket

read_arguments = False

if read_arguments:

    parser = argparse.ArgumentParser(
        description='Load saved joint and ball data from demonstrations. Process it and plot.')
    parser.add_argument('--joint_file', help='joint file')
    parser.add_argument('--ball_file', help='ball file')
    parser.add_argument(
        '--num_examples', help='number of demonstrations', type=int)
    parser.add_argument(
        '--plot_example', help='plot a specific example', type=int)
    parser.add_argument(
        '--smooth', help='smoothing factor of splines while plotting')
    parser.add_argument(
        '--align_with_ball', help='align racket with ball by a transformation')

    args = parser.parse_args()
    assert (args.plot_example <
            args.num_examples), "example to plot must be less than num of examples"
else:
    class MyArgs:
        joint_file = os.environ['HOME'] + \
            '/table-tennis/data/19.8.18/joints.txt'
        ball_file = os.environ['HOME'] + \
            '/table-tennis/data/19.8.18/balls.txt'
        num_examples = 4
        plot_example = 1
        smooth = 0.01
        align_with_ball = False

    args = MyArgs()

joint_file = args.joint_file  # './data/10.6.18/joints.txt'
ball_file = args.ball_file
num_examples = args.num_examples  # 21  # first dataset = 20, second = 21
M = np.fromfile(joint_file, sep=" ")
ndof = 7
dt = 0.002
N = M.size/(ndof+1)
M = np.reshape(M, [N, ndof+1])
q = M[:, 1:]
t_joints = M[:, 0]
t_joints = 0.001 * t_joints  # in miliseconds

t_min = t_joints[0]
if args.ball_file:
    B = np.fromfile(ball_file, sep=" ")
    N_balls = B.size/4
    B = np.reshape(B, [N_balls, 4])
    balls = B[:, 1:]
    # remove the zeros
    idx_nonzero = np.where(np.sum(balls, axis=1))[0]
    balls = balls[idx_nonzero, :]
    t_balls = B[idx_nonzero, 0]
    t_balls = 0.001 * t_balls
    t_min = min(t_min, t_balls[0])
    t_balls = t_balls - t_min
    t_joints = t_joints - t_min

'''
# quick plotting for ball
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(balls[:, 0], balls[:, 1], balls[:, 2], c="b")
plt.show()
'''


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
    # print clusters
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

    return low_vel_idxs.astype(int)


def plot_examples(examples, joint_dict, ball_dict=None, smooth_fact=0.01, align=False):

    idx_movements = joint_dict['idx_move']
    q = joint_dict['x']
    t = joint_dict['t']
    num_examples = examples.size
    for i in range(num_examples):
        f, axs = plt.subplots(7, 1, sharex=False)
        # print examples[i]
        idx_plot = np.arange(
            start=idx_movements[0, examples[i]],
            stop=idx_movements[1, examples[i]]+1, step=1, dtype=np.int32)
        # print idx_plot
        q_plot = q[idx_plot, :]
        t_plot = t[idx_plot]
        for j in range(7):
            spl = UnivariateSpline(
                t_plot, q_plot[:, j], w=None, k=3, s=smooth_fact)
            q_smooth = spl(t_plot)
            axs[j].plot(t_plot, q_plot[:, j])
            axs[j].plot(t_plot, q_smooth)
        # KINEMATICS PLOT
        if ball_dict is not None:
            x_plot = np.zeros((3, idx_plot.size))
            for idx, val in enumerate(idx_plot):
                As = wam.barrett_wam_kinematics(np.transpose(q_plot[idx, :]))
                R = As[-1, 0:3, 0:3]
                x_racket = As[-1, 0:3, -1]
                x_plot[:, idx] = x_racket

                if align:
                    quat = racket.rot2Quat(R)
                    orient = racket.calcRacketOrientation(quat)
                    R = np.squeeze(racket.quat2Rot(orient))
                    x_plot[:, idx] = x_plot[:, idx] + \
                        0.14 * R[:, 1] + 0.05 * R[:, 2]

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(x_plot[0, :], x_plot[1, :], x_plot[2, :], c="r")

            t_plot_robot = t_joints[idx_plot]
            idx_label_robot = np.arange(
                0, len(idx_plot), step=100, dtype=np.int32)
            for idx, x_robot in enumerate(x_plot[:, idx_label_robot]):
                label = str(t_plot_robot[idx_label_robot[idx]])
                ax.text(x_plot[0, idx_label_robot[idx]],
                        x_plot[1, idx_label_robot[idx]],
                        x_plot[2, idx_label_robot[idx]], label[:5])

            # extract ball
            balls = ball_dict['x']
            idx_move_ball = ball_dict['idx_move']
            t_balls = ball_dict['t']
            #t_balls = ball_dict['t']
            idx_plot = np.arange(
                idx_move_ball[0, examples[i]], idx_move_ball[1, examples[i]]+1, step=1, dtype=np.int32)
            balls_plot = balls[idx_plot, :]
            t_plot_ball = t_balls[idx_plot]
            #print balls_plot
            ax.scatter(balls_plot[:, 0], balls_plot[:, 1],
                       balls_plot[:, 2], c="b")
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            idx_label_ball = np.arange(
                0, len(idx_plot), step=10, dtype=np.int32)
            #print balls_plot[idx_label,:]

            for i, ball in enumerate(balls_plot[idx_label_ball, :]):
                label = str(t_plot_ball[idx_label_ball[i]])
                ax.text(balls_plot[idx_label_ball[i], 0],
                        balls_plot[idx_label_ball[i], 1],
                        balls_plot[idx_label_ball[i], 2], label[:5])
        plt.show()


idx_joint_move = detect_movements(q, num_examples, dt)
joint_dict = {'t': t_joints, 'x': q, 'idx_move': idx_joint_move}

# remove last
idx_joint_move = idx_joint_move[:, :-1]
num_examples = num_examples-1

# for each movement
# get the balls between t_joints
if args.ball_file:
    idx_ball_move = np.zeros((2, num_examples))
    for i in range(num_examples):
        idx_ball_move[0, i] = np.where(
            t_balls > t_joints[idx_joint_move[0, i]])[0][0]
        idx_ball_move[1, i] = np.where(
            t_balls < t_joints[idx_joint_move[1, i]])[0][-1]
    ball_dict = {'t': t_balls, 'x': balls, 'idx_move': idx_ball_move}
else:
    ball_dict = None

if args.smooth:
    smooth = args.smooth
else:
    smooth = 0.01

# examples = np.array([8])
if args.plot_example is not None:
    examples = np.array([args.plot_example])
    plot_examples(examples, joint_dict, ball_dict,
                  smooth_fact=smooth, align=args.align_with_ball)
