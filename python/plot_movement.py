'''
load data
cut it and process
plot an example
'''

import numpy as np

filename = './data/10.6.18/joints.txt'
M = np.fromfile(filename, sep=" ")
ndof = 7
N = M.size/ndof
q = np.reshape(M, [N, ndof])
dt = 0.002
t = dt * np.linspace(1, N, N)
num_examples = 21  # first dataset = 20, second = 21


def detect_movements(q, num_examples, dt):

    # Find the prominent fast moving segments
    qd = np.diff(q, 1, 0)/dt
    vels = np.sqrt(np.sum(qd*qd, -1))
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


idx_movements = detect_movements(q, num_examples, dt)
#print idx_movements


def plot_examples(examples, t, q, idx_movements):
    import matplotlib.pyplot as plt
    from scipy.interpolate import UnivariateSpline
    #from mpl_toolkits.mplot3d import Axes3D

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
        plt.show()


# plot one of the motions
'''
GOOD MOTIONS (21 TOTAL)
[2,4,6,9,10,12,15,16,18,19,20,21]-1 (zero indexing)
'''
examples = np.array([8])
plot_examples(examples, t, q, idx_movements)
