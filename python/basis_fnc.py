'''
Gaussian Basis functions and their all sorts of derivatives
are stored here
'''

import numpy as np


def create_acc_weight_mat(centers, widths, time, include_intercept=True):
    '''
    Create weighting matrix that penalizes the accelerations
    '''
    num_bumps = len(centers)
    N = len(time)
    if include_intercept:
        Xdot2 = np.zeros((N, num_bumps+1))  # also include intercept
    else:
        Xdot2 = np.zeros((N, num_bumps))  # also include intercept
    for i in range(num_bumps):
        Xdot2[:, i] = basis_fnc_gauss_der2(time, centers[i], widths[i])

    return np.dot(Xdot2.T, Xdot2), Xdot2


def create_acc_weight_mat_der(centers, widths, time, include_intercept=True, der='c'):

    num_bumps = len(centers)
    N = len(time)
    if include_intercept:
        Xdot2 = np.zeros((N, num_bumps+1))  # also include intercept
    else:
        Xdot2 = np.zeros((N, num_bumps))  # also include intercept

    for i in range(num_bumps):
        if der == 'c':
            Xdot2[:, i] = basis_fnc_gauss_der2_der_c(
                time, centers[i], widths[i])
        elif der == 'w':
            Xdot2[:, i] = basis_fnc_gauss_der2_der_w(
                time, centers[i], widths[i])
        else:
            raise Exception('Derivative not recognized!')
    return Xdot2


def create_gauss_regressor(centers, widths, time, include_intercept=True):
    '''
    Create Gaussian basis functions num_bumps many times along
    '''
    num_bumps = len(centers)
    N = len(time)
    if include_intercept:
        X = np.zeros((N, num_bumps+1))  # also include intercept
    else:
        X = np.zeros((N, num_bumps))

    for i in range(num_bumps):
        X[:, i] = basis_fnc_gauss(time, centers[i], widths[i])

    if include_intercept:
        X[:, -1] = 1.0
    return X


def create_gauss_regressor_der(centers, widths, time, include_intercept=True, der='c'):
    num_bumps = len(centers)
    N = len(time)
    if include_intercept:
        X = np.zeros((N, num_bumps+1))  # also include intercept
    else:
        X = np.zeros((N, num_bumps))

    if der == 'w':
        for i in range(num_bumps):
            X[:, i] = basis_fnc_gauss_der_w(time, centers[i], widths[i])
    elif der == 'c':
        for i in range(num_bumps):
            X[:, i] = basis_fnc_gauss_der_c(time, centers[i], widths[i])
    else:
        raise Exception('derivative not known!')

    if include_intercept:
        X[:, -1] = 0.0
    return X


def basis_fnc_gauss(t, c, w):

    return np.exp(-((t-c)**2)/w)


def basis_fnc_gauss_der2(t, c, w):
    ''' 2nd derivative w.r.t time '''
    return (-(2.0/w) + ((4.0/(w*w))*(t-c)**2))*basis_fnc_gauss(t, c, w)


def basis_fnc_gauss_der_c(t, c, w):
    ''' 1st der w.r.t c '''
    return 2*(t-c)/w * basis_fnc_gauss(t, c, w)


def basis_fnc_gauss_der_w(t, c, w):
    ''' 1st der w.r.t widths w'''
    return ((t-c)/w)**2 * basis_fnc_gauss(t, c, w)


def basis_fnc_gauss_der2_der_c(t, c, w):
    '''c-derivative of 2nd time derivative of basis fnc'''
    return -(4.0/(w*w))*2*(t-c)*basis_fnc_gauss(t, c, w) + \
        (-2/w + (4.0/(w*w))*(t-c)**2) * basis_fnc_gauss_der_c(t, c, w)


def basis_fnc_gauss_der2_der_w(t, c, w):
    '''w-derivative of 2nd time der. of basis fnc'''
    return (2.0/(w*w) - ((8/(w*w*w))*(t-c)**2))*basis_fnc_gauss(t, c, w) + \
        (-(2.0/w) + ((4.0/(w*w))*(t-c)**2))*basis_fnc_gauss_der_w(t, c, w)
