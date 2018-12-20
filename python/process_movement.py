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
#from matplotlib2tikz import save as tikz_save


def load_files(args):

    joint_file = args.joint_file  # './data/10.6.18/joints.txt'
    ball_file = args.ball_file
    M = np.loadtxt(joint_file)
    ndof = 7
    if M.shape[1] == (ndof + 1):
        q = M[:, 1:]
        t_joints = M[:, 0]
        t_joints = 0.001 * t_joints  # in miliseconds
    elif M.shape[1] == ndof:  # older dataset
        N = M.shape[0]
        q = M  # np.reshape(M, [N, ndof])
        t_joints = 0.002 * np.linspace(0, N, N)
    else:
        raise Exception('File doesnt have right size!')

    t_balls = []
    balls = []
    if args.ball_file is not None:
        B = np.fromfile(ball_file, sep=" ")
        N_balls = B.size/4
        B = np.reshape(B, [N_balls, 4])
        balls = B[:, 1:]
        # remove the zeros
        idx_nonzero = np.where(np.sum(balls, axis=1))[0]
        balls = balls[idx_nonzero, :]
        t_balls = B[idx_nonzero, 0]
        t_balls = 0.001 * t_balls
        t_min = min(t_joints[0], t_balls[0])
        t_balls = t_balls - t_min
        t_joints = t_joints - t_min

    return t_joints, t_balls, q, balls


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


def compute_kinematics(idx, q, align):
    '''
    Compute kinematics (x positions) of racket center positions
    Align racket center with ball by transforming 14 cm back (y) and 5 cm in z-dir
    '''
    x_plot = np.zeros((3, idx.size))
    for idx, val in enumerate(idx):
        As = wam.barrett_wam_kinematics(np.transpose(q[idx, :]))
        R = As[-1, 0:3, 0:3]
        x_racket = As[-1, 0:3, -1]
        x_plot[:, idx] = x_racket

        if align:
            quat = racket.rot2Quat(R)
            orient = racket.calcRacketOrientation(quat)
            R = np.squeeze(racket.quat2Rot(orient))
            x_plot[:, idx] = x_plot[:, idx] + \
                0.14 * R[:, 1] + 0.05 * R[:, 2]
    return x_plot


def plot_joints(t_plot, q_plot, smooth_opts=None):
    f, axs = plt.subplots(7, 1, sharex=False)
    for j in range(7):
        axs[j].plot(t_plot, q_plot[:, j])
        if smooth_opts is not None:
            w = smooth_opts.weights
            s = smooth_opts.factor
            k = smooth_opts.degree
            spl = UnivariateSpline(
                t_plot, q_plot[:, j], w=w, k=k, s=s)
            q_smooth = spl(t_plot)
            axs[j].plot(t_plot, q_smooth)


def plot_3d_with_ball(idx_plot, q_plot, t_joints, ball_dict, align=False):
    x_plot = compute_kinematics(idx_plot, q_plot, align)
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
    t_balls = ball_dict['t']
    # print balls_plot
    ax.scatter(balls[:, 0], balls[:, 1], balls[:, 2], c="b")
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    idx_label_ball = np.arange(
        0, balls.shape[0], step=10, dtype=np.int32)
    # print balls_plot[idx_label,:]

    for i, ball in enumerate(balls[idx_label_ball, :]):
        label = str(t_balls[idx_label_ball[i]])
        ax.text(balls[idx_label_ball[i], 0],
                balls[idx_label_ball[i], 1],
                balls[idx_label_ball[i], 2], label[:5])
    return x_plot


def plot_example(example, joint_dict, ball_dict=None, smooth_opts=None, align=False, dump=True):
    ''' Plots an example demonstration in joint space and in task space '''

    idx_movements = joint_dict['idx_move']
    q = joint_dict['x']
    t_joints = joint_dict['t']
    # print examples[i]
    idx_plot = np.arange(
        start=idx_movements[0, example],
        stop=idx_movements[1, example]+1, step=1, dtype=np.int32)
    # print idx_plot
    q_plot = q[idx_plot, :]
    t_plot = t_joints[idx_plot]
    plot_joints(t_plot, q_plot, smooth_opts)
    # KINEMATICS PLOT
    if ball_dict is not None:
        x_plot = plot_3d_with_ball(
            idx_plot, q_plot, t_joints, ball_dict, align=align)
    plt.show(block=False)

    if dump:
        # dump to a text file
        robot_mat = np.hstack((t_plot[:, np.newaxis], q_plot))
        np.savetxt('python/dumps/move_robot_' + str(example) + '.txt',
                   X=robot_mat, fmt='%.4e', delimiter='\t')
        if ball_dict is not None:
            robot_mat = np.hstack((robot_mat, x_plot.T))
            ball_mat = np.hstack(
                (ball_dict['t'][:, np.newaxis], ball_dict['x']))
            np.savetxt('python/dumps/move_ball_' + str(example) + '.txt',
                       X=ball_mat, fmt='%.4e', delimiter='\t')


def get_ball_obs_for_movement(t_balls, balls, t_joints, idx_joint_move, remove_outlier=True):

    idx_ball_move = np.zeros((2,), dtype=np.int64)
    idx_ball_move[0] = int(np.where(
        t_balls > t_joints[idx_joint_move[0]])[0][0])
    idx_ball_move[1] = int(np.where(
        t_balls < t_joints[idx_joint_move[1]])[0][-1])

    ''' if remove_outlier is turned on, then removes balls that jump'''
    # print idx_ball_move

    balls = balls[idx_ball_move[0]:idx_ball_move[1], :]
    if remove_outlier:
        len_idx = idx_ball_move[1]-idx_ball_move[0]
        idx_robust = np.zeros((len_idx,), dtype=np.int64)
        idx_robust[0] = 0
        thresh = 0.1  # 10 cm diff.
        k = 1
        for i in range(len_idx-1):
            if np.linalg.norm(balls[i+1, :]-balls[i, :]) < thresh:
                idx_robust[k] = i+1
                k += 1
        balls = balls[idx_robust, :]
        t_balls = t_balls[idx_robust]

    ball_dict = {'t': t_balls, 'x': balls, 'idx_move': idx_ball_move}
    return ball_dict


def run_serve_demo(args):
    ''' Process one serve demonstration. Main entrance point to process_movement.

    Also plots the data if there is a demand.
    '''

    t_joints, t_balls, q, balls = load_files(args)
    idx_joint_move = detect_movements(q, args.num_examples, dt=0.002)
    joint_dict = {'t': t_joints, 'x': q, 'idx_move': idx_joint_move}

    # for each movement
    # get the balls between t_joints
    if args.ball_file:
        ball_dict = get_ball_obs_for_movement(
            t_balls, balls, t_joints, idx_joint_move[:, args.plot_example])
    else:
        ball_dict = None

    # examples = np.array([8])
    if args.plot:
        plot_example(args.plot_example, joint_dict, ball_dict,
                     smooth_opts=args.smooth, align=args.align_with_ball, dump=args.dump_plot_data)
    return joint_dict, ball_dict


def process_args():
    '''
    @deprecated! Process input arguments
    '''

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
    return args


def create_default_args(date='24.8.18'):

    class MyArgs:
        joint_file = os.environ['HOME'] + \
            '/table-tennis/data/' + date + '/joints.txt'
        ball_file = os.environ['HOME'] + \
            '/table-tennis/data/' + date + '/balls.txt'
        num_examples = 4
        plot_example = 0  # TODO: ball_dict returns only plot_example I think
        plot = True
        dump_plot_data = False

        class Smooth:
            factor = 0.01
            weights = None
            degree = 3
        smooth = Smooth()
        align_with_ball = False
    return MyArgs()


def save_to_tikz(path='../learning-to-serve/Pictures/', name='demo_kin.tex'):
    # axis_opts = {'xlabel_near_ticks', 'ylabel_near_ticks', 'scale_only_axis'}
    tikz_save(path+name, figureheight='\\figureheight',
              figurewidth='\\figurewidth')  # extra_axis_parameters=axis_opts)


if __name__ == '__main__':
    # args = process_movement()
    # date = '15.11.18'
    # args = create_default_args(date)
    # args.num_examples = 15
    # args.plot_example = 5
    # args.ball_file = None
    # args.smooth = None  # dont draw smoothened version
    # args.dump_plot_data = False  # set to true to dump plot data
    args = create_default_args()
    run_serve_demo(args)
