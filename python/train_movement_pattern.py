'''
Train movement pattern from examples

'''
import sys
sys.path.append('python/')
import numpy as np
import process_movement as serve
from sklearn import linear_model as lm
import matplotlib.pyplot as plt
import json


def plot_regression(X, theta, q):
    q_smooth = np.dot(X, theta)
    f, axs = plt.subplots(7, 1, sharex=True)
    t = 0.002 * np.arange(1, X.shape[0]+1)
    for i in range(7):
        axs[i].plot(t, q[:, i])
        axs[i].plot(t, q_smooth[:, i])
    plt.show()


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


def basis_fnc_gauss(t, c, w):

    return np.exp(-((t-c)**2)/w)


def basis_fnc_gauss_der2(t, c, w):

    return (-(2.0/w) + ((4.0/(w*w))*(t-c)**2))*basis_fnc_gauss(t, c, w)


def train_sparse_weights(plot_regr_results=False):

    args = serve.create_default_args()  # load default options
    args.plot = False
    args.proces_example = 0  # change example
    '''
    args.smooth.factor = 0.01
    args.smooth.degree = 3
    args.smooth.extrap = 3  # doesnt seem to change anything
    '''
    joint_dict, ball_dict = serve.run_serve_demo(args)

    # train multi task elastic net
    example = 0
    idx_move = joint_dict['idx_move']
    idx = np.arange(idx_move[0, example], idx_move[1, example])
    q = joint_dict['x'][idx, :]
    t = joint_dict['t'][idx]
    p = 100  # number of features in total
    c = np.linspace(0, t[-1]-t[0], p)  # centers
    w = 0.05 * t[-1] * np.ones((p,))  # widths
    X = create_gauss_regressor(c, w, time=t-t[0])

    # theta = np.linalg.lstsq(X, q)[0]
    # clf = lm.Lasso(alpha=0.001, fit_intercept=False)
    # clf = lm.MultiTaskLasso(alpha=0.0024, fit_intercept=False)
    # running multi task laso with cross-validation gives 0.002
    # clf = lm.MultiTaskLassoCV(eps=1e-3, n_alphas=100, alphas=None,
    #                                    fit_intercept=False, cv=10, verbose=True, n_jobs=-1)

    # l1_ratio_list = np.linspace(0.1, 1.0, 10)
    # l1_ratio_list = 1-np.exp(-np.arange(1, 10)/2.0)
    # clf = lm.MultiTaskElasticNetCV(l1_ratio=l1_ratio_list, eps=1e-3, n_alphas=100, alphas=None,
    #                               fit_intercept=False, cv=3, verbose=True, n_jobs=-1)

    # running cross-val gives alpha = 0.0038, l1_ratio = 0.632
    # clf = lm.MultiTaskElasticNet(
    #     alpha=0.0038, l1_ratio=0.632, fit_intercept=False)
    # clf.fit(X, q)

    # transform multitask elastic net to multitask lasso

    alpha = 0.002
    rho = 0.99
    N = q.shape[0]
    d = q.shape[1]  # 7
    lamb1 = 2*N*alpha*rho
    lamb2 = N*alpha*(1-rho)
    q_bar = np.vstack((q, np.zeros((N, d))))
    mult = np.sqrt(1.0/(1+lamb2))
    _, Xdot2 = create_acc_weight_mat(c, w, t, include_intercept=True)
    # print Xdot2[0:5, 0:5]
    X_bar = mult * np.vstack((X, np.sqrt(lamb2)*Xdot2))
    # clf = lm.MultiTaskLassoCV(eps=1e-3, n_alphas=100, alphas=None,
    #                          fit_intercept=False, cv=10, verbose=True, n_jobs=-1)
    clf = lm.MultiTaskLasso(alpha=alpha, fit_intercept=False)
    clf.fit(X_bar, q_bar)
    theta = clf.coef_.T * mult
    # theta = clf.coef_.T

    if plot_regr_results:
        plot_regression(X, theta, q)

    idx_non = np.nonzero(theta[:, 0])[0]  # only nonzero entries
    theta = theta[idx_non, :]
    # last one params are the intercepts!
    dump_json_regression_obj(c[idx_non[:-1]], w[idx_non[:-1]], theta)
    # return X[:, idx_non], theta[idx_non, :]


def dump_json_regression_obj(centers, widths, theta, basis_type="squared exp"):
    '''
    Create a dictionary and dump it as json object.
    '''
    json_regr = dict()
    json_regr['basis_type'] = basis_type
    json_regr['intercepts'] = theta[-1, :].tolist()
    json_regr['centers'] = centers.tolist()
    json_regr['widths'] = widths.tolist()
    joint_params = list()
    for i in range(7):
        params = dict()
        params['ID'] = i
        params['params'] = theta[:-1, i].tolist()
        joint_params.append(params)
    json_regr['joints'] = joint_params
    file_name = 'json/rbf.json'
    with open(file_name, "w") as write_file:
        json.dump(json_regr, write_file)
