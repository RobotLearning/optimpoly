'''
load data
cut it and process
plot an example
'''

import numpy as np
import argparse
import barrett_wam_kinematics as wam

parser = argparse.ArgumentParser(
    description='Load saved joint and ball data from demonstrations. Process it and plot.')
parser.add_argument('--joint_file', help='joint file')
parser.add_argument('--ball_file', help='ball file')
parser.add_argument(
    '--num_examples', help='number of demonstrations', type=int)
parser.add_argument('--plot_example', help='plot a specific example', type=int)
parser.add_argument(
    '--plot_3d', help='also plot the cartesian values', action="store_true")
args = parser.parse_args()
assert (args.plot_example <
        args.num_examples), "example to plot must be less than num of examples"
joint_file = args.joint_file  # './data/10.6.18/joints.txt'
num_examples = args.num_examples  # 21  # first dataset = 20, second = 21
M = np.fromfile(joint_file, sep=" ")
ndof = 7
N = M.size/ndof
q = np.reshape(M, [N, ndof])
dt = 0.002
t = dt * np.linspace(1, N, N)

if args.ball_file:
    B = np.fromfile(ball_file, sep=" ")
    N_balls = B.size/3
    balls = np.reshape(B, [N_balls, 3])
    t_balls = dt * np.linspace(1, N, N)


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
    thresh = 1000  # 2 sec diff min betw demonstr.
    while idx_clusters.size < num_examples:
        diff_max = min(abs(idx_vels[idx] - idx_clusters))
        if diff_max > thresh:
            idx_clusters = np.insert(idx_clusters, 0, idx_vels[idx])
        idx = idx+1

    # sort the points and find points of low velocity around them
    clusters = np.sort(idx_clusters)
    low_vel_idxs = np.zeros((2, num_examples))
    low_vel_thresh = 1e-3
    max_duration_movement = 1.0  # seconds
    idx_max_move = max_duration_movement/dt

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


def plot_examples(examples, t, q, idx_movements, plot_3d=False):
    import matplotlib.pyplot as plt
    from scipy.interpolate import UnivariateSpline
    from mpl_toolkits.mplot3d import Axes3D

    num_examples = examples.size
    for i in range(num_examples):
        f, axs = plt.subplots(7, 1, sharex=False)
        # print examples[i]
        idx_plot = np.arange(
            start=idx_movements[0, examples[i]],
            stop=idx_movements[1, examples[i]]+1, step=1, dtype=np.int32)
        q_plot = q[idx_plot, :]
        t_plot = t[idx_plot]
        for j in range(7):
            spl = UnivariateSpline(
                t_plot, q_plot[:, j], w=None, k=3, s=0.2)
            q_smooth = spl(t_plot)
            axs[j].plot(t_plot, q_plot[:, j])
            axs[j].plot(t_plot, q_smooth)
        # KINEMATICS PLOT
        if plot_3d:
            x_plot = np.zeros((3, idx_plot.size))
            for idx, val in enumerate(idx_plot):
                As = wam.barrett_wam_kinematics(np.transpose(q_plot[idx, :]))
                x_plot[:, idx] = np.transpose(As[-1])[-1][0:3]
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            for i in range(3):
                ax.scatter(x_plot[0, :], x_plot[1, :], x_plot[2, :], c="b")
        plt.show()


# plot one of the motions
'''
GOOD MOTIONS (21 TOTAL)
[2,4,6,9,10,12,15,16,18,19,20,21]-1 (zero indexing)
'''
idx_joint_movements = detect_movements(q, num_examples, dt)
idx_ball_movements = detect_movements(balls, num_examples, dt)
#print idx_movements


#examples = np.array([8])
if args.plot_example:
    examples = np.array([args.plot_example])
    if args.plot_3d:
        plot_examples(examples, t, q, idx_joint_movements, plot_3d=True)
        if args.ball_file:
            plot_examples(examples, t_balls, balls,
                          idx_ball_movements, plot_3d=True)
    else:
        plot_examples(examples, t, q, idx_joint_movements)
