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
from scipy import optimize as opt


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


def dump_json_regression_obj(centers, widths, theta, basis_type="squared exp", file_name="rbf.json"):
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
    file_name = 'json/' + file_name
    print 'Saving to ', file_name
    with open(file_name, "w") as write_file:
        json.dump(json_regr, write_file)


def l2_pen_regr(X, q, C):
    ''' L2 penalized standard regression'''
    # theta = np.linalg.lstsq(X, q)[0]
    M = np.dot(X.T, X) + 1e-5*C
    theta = np.linalg.solve(M, np.dot(X.T, q))
    res = q - np.dot(X, theta)
    return theta, res


def multi_task_lasso(X, q, cv=False, alpha=0.002):
    '''Multi Task Lasso with dimensions forced to share features
    Running multi task Lasso with cross-validation gives 0.002
    '''
    if cv:
        clf = lm.MultiTaskLassoCV(eps=1e-3, n_alphas=100, alphas=None,
                                  fit_intercept=False, cv=10, verbose=True, n_jobs=-1)
    else:
        clf = lm.MultiTaskLasso(alpha=alpha, fit_intercept=False)
    clf.fit(X, q)
    theta = clf.coef_.T
    res = q - np.dot(X, theta)
    return theta, res


def multi_task_elastic_net(X, q, cv=False, alpha=0.0038, l1_ratio=0.632):
    '''
    Multi Task Elastic Net with dimensions forced to share features
    both l1 and l2 regularization is employed in the Elastic Net formulation

    Running cross-val gives alpha = 0.0038, l1_ratio = 0.632
    '''
    if cv:
        l1_ratio_list = np.linspace(0.1, 1.0, 10)
        l1_ratio_list = 1-np.exp(-np.arange(1, 10)/2.0)
        clf = lm.MultiTaskElasticNetCV(l1_ratio=l1_ratio_list, eps=1e-3, n_alphas=100, alphas=None,
                                       fit_intercept=False, cv=3, verbose=True, n_jobs=-1)
    else:
        clf = lm.MultiTaskElasticNet(
            alpha=alpha, l1_ratio=l1_ratio, fit_intercept=False)
    clf.fit(X, q)
    theta = clf.coef_.T
    res = q - np.dot(X, theta)
    return theta, res


def multi_task_weighted_elastic_net(X, q, L, alpha=0.001, rho=0.99999):
    '''
    Since sklearn does not support weighted version (for us weighted L2 regularization) of Elastic Net, we transform multitask elastic net to multitask lasso
    solve it and then transform the solution back to elastic net

    L is the accelerations that composes the weight matrix, i.e. W = L'*L
    '''
    N = q.shape[0]
    d = q.shape[1]  # 7
    lamb1 = 2*N*alpha*rho
    lamb2 = N*alpha*(1-rho)
    q_bar = np.vstack((q, np.zeros((N, d))))
    mult = np.sqrt(1.0/(1+lamb2))
    Xdot2 = L  # create_acc_weight_mat(c, w, t, include_intercept=True)
    # print Xdot2[0:5, 0:5]
    X_bar = mult * np.vstack((X, np.sqrt(lamb2)*Xdot2))
    # clf = lm.MultiTaskLassoCV(eps=1e-3, n_alphas=100, alphas=None,
    #                          fit_intercept=False, cv=10, verbose=True, n_jobs=-1)
    clf = lm.MultiTaskLasso(alpha=alpha, fit_intercept=False)
    clf.fit(X_bar, q_bar)
    theta = clf.coef_.T * mult
    res = q - np.dot(X, theta)

    return theta, res


def train_sparse_weights(plot_regr_results=False, ex=0, save=False):

    args = serve.create_default_args()  # load default options
    args.plot = False
    args.process_example = ex  # change example
    '''
    args.smooth.factor = 0.01
    args.smooth.degree = 3
    args.smooth.extrap = 3  # doesnt seem to change anything
    '''
    joint_dict, ball_dict = serve.run_serve_demo(args)

    # train multi task elastic net
    example = args.process_example
    idx_move = joint_dict['idx_move']
    idx = np.arange(idx_move[0, example], idx_move[1, example])
    q = joint_dict['x'][idx, :]
    t = joint_dict['t'][idx]
    t -= t[0]
    p = 100  # number of features in total
    c = np.linspace(t[0], t[-1], p)  # centers
    w = 0.1 * np.ones((p,))  # widths
    X = create_gauss_regressor(c, w, time=t)
    C, Xdot2 = create_acc_weight_mat(c, w, t, include_intercept=True)  #

    # theta, res = l2_pen_regr(X, q, C)
    # theta, res = multi_task_lasso(X, q)
    # theta, res = multi_task_elastic_net(X, q)
    # theta, res = multi_task_weighted_elastic_net(X, q, Xdot2)
    X, theta, c, w = iter_lasso(t, q)
    if plot_regr_results:
        plot_regression(X, theta, q)

    idx_non = np.nonzero(theta[:, 0])[0]  # only nonzero entries
    print 'Lasso kept', len(idx_non), 'solutions!'
    theta = theta[idx_non, :]
    # last one params are the intercepts!
    if save:
        filename = "rbf_" + str(example) + ".json"
        dump_json_regression_obj(
            c[idx_non[:-1]], w[idx_non[:-1]], theta, file_name=filename)
    # return X[:, idx_non], theta[idx_non, :]


def iter_lasso(t, q):

    # initialize the iteration
    p = 400
    c = np.linspace(t[0], t[-1], p) + 0.01 * np.random.randn(p)
    w = 0.1 * np.ones((p,)) + 0.01 * np.random.rand(p)
    X = create_gauss_regressor(c, w, t)
    _, Xdot2 = create_acc_weight_mat(c, w, t)
    iter_max = 3
    a = 0.001  # alpha
    r = 0.99999  # rho
    N = q.shape[0]
    lamb1 = 2*N*a*r
    lamb2 = N*a*(1-r)

    def f_opt(x):
        c_opt = x[:p]
        w_opt = x[p:]
        X_new = create_gauss_regressor(c_opt, w_opt, t)
        C, _ = create_acc_weight_mat(c_opt, w_opt, t)
        res = q - np.dot(X_new, theta)
        cost = np.linalg.norm(res, 'fro')
        theta_21_norm = np.sum(np.linalg.norm(theta, axis=1))
        l2_acc_pen = np.linalg.norm(np.dot(theta.T, np.dot(C, theta)), 'fro')
        cost += lamb1 * theta_21_norm
        cost += lamb2 * l2_acc_pen
        return cost

    def prune_params(params):  # acts globally in the function
        idx_non = np.nonzero(params[:, 0])[0]  # only nonzero entries
        p_new = len(idx_non)-1
        print 'Lasso kept', p_new, 'solutions!'  # intercept is static
        theta_new = params[idx_non, :]
        c_new = c[idx_non[:-1]]
        w_new = w[idx_non[:-1]]
        return theta_new, c_new, w_new, p_new

    xopt = np.vstack((c, w))
    theta, residual = multi_task_weighted_elastic_net(
        X, q, Xdot2, alpha=a, rho=r)

    for i in range(iter_max):
        theta, c, w, p = prune_params(theta)
        # update RBF weights
        result = opt.minimize(f_opt, xopt, method="BFGS")
        xopt = result.x
        c = xopt[:p]
        w = xopt[p:]
        # print c, w
        X = create_gauss_regressor(c, w, t)
        _, Xdot2 = create_acc_weight_mat(c, w, t)
        # perform lasso
        res_last = residual
        theta, residual = multi_task_weighted_elastic_net(X, q, Xdot2)
        # shrink the regularizers
        # to scale lasso throughout iterations
        a /= np.linalg.norm(res_last)/np.linalg.norm(residual)
        lamb1 = 2*N*a*r
        lamb2 = N*a*(1-r)

    return X, theta, c, w


if __name__ == '__main__':
    train_sparse_weights(True, ex=2, save=True)
