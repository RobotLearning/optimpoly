import sys
import time
sys.path.append('python/')
import numpy as np
import process_movement as serve
import basis_fnc as basis
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
    plt.show(block=False)


def plot_multi_regression(X, theta, Q, intercepts):
    ''' X is (7N x M) matrix, theta is (MxD), Q is (7N x D) '''
    Q_smooth = np.dot(X, theta)

    N = X.shape[0]/7
    D = theta.shape[1]
    t = 0.002 * np.arange(1, N+1)
    for d in range(D):
        q = Q[:, d]
        q = np.reshape(q, ((7, N)))
        q_smooth = Q_smooth[:, d]
        q_smooth = np.reshape(q_smooth, ((7, N))) + \
            intercepts[:, d][:, np.newaxis]
        f, axs = plt.subplots(7, 1, sharex=True)
        for i in range(7):
            axs[i].plot(t, q[i, :])
            axs[i].plot(t, q_smooth[i, :])
    plt.show(block=False)


def dump_json_regression_obj(centers, widths, theta, basis_type="squared exp", file_name="rbf.json", order=None):
    '''
    Create a dictionary and dump it as json object.
    '''
    json_regr = dict()
    json_regr['basis_type'] = basis_type
    json_regr['centers'] = centers.tolist()
    json_regr['widths'] = widths.tolist()
    joint_params = list()
    for i in range(7):
        params = dict()
        params['ID'] = i
        params['intercept'] = theta[-1, i]
        params['params'] = theta[:-1, i].tolist()
        joint_params.append(params)
    json_regr['joints'] = joint_params
    if order is not None:
        json_regr['order'] = order
    file_name = 'json/' + file_name
    print 'Saving to ', file_name
    with open(file_name, "w") as write_file:
        json.dump(json_regr, write_file)


def dump_json_multi_demo_regression_obj(centers, widths, theta, intercepts, examples, args, basis_type="squared exp", order=None):
    '''
    Create a dictionary and dump it as json object.
    '''
    ndof = 7
    date = args.date
    num_basis_per_joint = len(centers)/ndof
    for j in range(len(examples)):
        json_regr = dict()
        json_regr['basis_type'] = basis_type
        json_regr['params'] = theta[:, j].tolist()
        joint_params = list()
        for i in range(7):
            idx = np.arange(i*num_basis_per_joint, (i+1)
                            * num_basis_per_joint, 1)
            params = dict()
            params['ID'] = i
            params['centers'] = centers[idx].tolist()
            params['widths'] = widths[idx].tolist()
            params['intercept'] = intercepts[i, j]
            joint_params.append(params)
        json_regr['joints'] = joint_params
        if order is not None:
            json_regr['order'] = order
        file_name = 'json/rbf_' + str(examples[j]) + '_' + date + '.json'
        print 'Saving to ', file_name
        with open(file_name, "w") as write_file:
            json.dump(json_regr, write_file)


def train_l2_reg_regr(args, plot_regr=True, ex=18, save=False, p=10, verbose=True):
    ''' Here we have multiple demonstrations '''
    joint_dict, ball_dict = serve.run_serve_demo(args)
    # train multi task elastic net
    idx_move = joint_dict['idx_move']
    idx = np.arange(idx_move[0, ex], idx_move[1, ex])
    q = joint_dict['x'][idx, :]
    t = joint_dict['t'][idx]
    t -= t[0]
    c = np.linspace(t[0], t[-1], p)  # centers
    w = 0.1 * np.ones((p,))  # widths
    X = basis.create_gauss_regressor(c, w, time=t)
    C, Xdot2 = basis.create_acc_weight_mat(c, w, t, include_intercept=True)  #

    theta, res = l2_pen_regr(t, X, q, C, 1e-4)
    if verbose:
        print 'Res. norm:', np.linalg.norm(q - np.dot(X, theta), 'fro')
        print 'Acc. norm:', np.linalg.norm(np.dot(Xdot2, theta), 'fro')
        print 'No. params:', theta.size


def l2_pen_regr(t, X, q, C, lamb):
    ''' L2 penalized standard regression'''
    # theta = np.linalg.lstsq(X, q)[0]
    M = np.dot(X.T, X) + lamb*C
    theta = np.linalg.solve(M, np.dot(X.T, q))
    res = q - np.dot(X, theta)
    return theta, res


def iter_multi_dof_lasso(t, q, p):
    '''
    Multi-Task grouping is performed across degrees of freedom(joints).
    Iterative MultiTaskElasticNet with nonlinear optimization(BFGS) to update BOTH the RBF
    parameters as well as the regression parameters.
    '''
    import multi_dof_lasso as lasso
    # initialize the iteration
    c = np.linspace(t[0], t[-1], p) + 0.01 * np.random.randn(p)
    w = 0.1 * np.ones((p,)) + 0.01 * np.random.rand(p)
    X = basis.create_gauss_regressor(c, w, t)
    _, Xdot2 = basis.create_acc_weight_mat(c, w, t)
    iter_max = 3
    a = 0.001  # alpha
    r = 0.99999  # rho
    N = q.shape[0]
    lamb1 = 2*N*a*r
    lamb2 = N*a*(1-r)

    def f_opt(x):
        c_opt = x[:p]
        w_opt = x[p:]
        f = lasso.elastic_net_cost(c_opt, w_opt, t, q, theta, lamb1, lamb2)
        df = lasso.elastic_net_cost_der(c_opt, w_opt, t, q, theta, lamb2)
        return f, df

    xopt = np.hstack((c, w))
    theta, residual, _ = lasso.multi_task_weighted_elastic_net(
        X, q, Xdot2, alpha=a, rho=r)

    for i in range(iter_max):
        theta, c, w, p = lasso.prune_params(theta, c, w)
        xopt = np.hstack((c, w))
        # update RBF weights
        result = opt.minimize(f_opt, xopt, jac=True, method="BFGS")
        xopt = result.x
        c = xopt[:p]
        w = xopt[p:]
        # print c, w
        X = basis.create_gauss_regressor(c, w, t)
        _, Xdot2 = basis.create_acc_weight_mat(c, w, t)
        # perform lasso
        res_last = residual
        theta, residual, _ = lasso.multi_task_weighted_elastic_net(
            X, q, Xdot2, alpha=a, rho=r)
        # shrink the regularizers
        # to scale lasso throughout iterations
        a /= (np.linalg.norm(res_last, 'fro') /
              np.linalg.norm(residual, 'fro'))**2
        lamb1 = 2*N*a*r
        lamb2 = N*a*(1-r)

    theta, c, w, p = lasso.prune_params(theta, c, w)
    X = basis.create_gauss_regressor(c, w, t)
    Xdot2 = basis.create_gauss_regressor_der(c, w, t)
    return X, theta, c, w, a, r, Xdot2


def iter_multi_demo_lasso(t, Q, p, include_intercept=False):
    '''
    Multi-Task grouping is performed across multiple demonstrations!
    Iterative MultiTaskElasticNet with nonlinear optimization(BFGS) to update BOTH the RBF
    parameters as well as the regression parameters.
    In the case of Barrett WAM there are 7x more RBF parameters than
    iter_multi_dof_lasso!

    Q is a IxJxK 3d-array
    t is a I-vector
    where I =  # of time points
          J =  # of dofs
          K =  # of demos

    An intercept is also learned, X will have a column of 1's
    '''
    import multi_dof_lasso as lasso
    import multi_demo_lasso as lasso_demo
    reload(lasso_demo)

    # initialize the iteration
    n_dofs = 7  # degrees of freedom
    n_tp = Q.shape[0]/n_dofs  # number of time points
    n_demos = Q.shape[1]

    # initialize big X
    c = np.linspace(t[0], t[-1], p) + 0.01 * np.random.randn(p)
    w = 0.1 * np.ones((p,)) + 0.01 * np.random.rand(p)
    c = np.tile(c, (n_dofs,))
    w = np.tile(w, (n_dofs,))
    X, Xdot2 = lasso_demo.stack_regressors(c, w, t, n_dofs, include_intercept)
    iter_max = 3
    a = 0.0001  # alpha
    r = 0.99999  # rho
    N = n_tp * n_dofs
    lamb1 = 2*N*a*r
    lamb2 = N*a*(1-r)

    def f_opt(x):
        c_new = x[:n_dofs*p]
        w_new = x[n_dofs*p:]
        f = lasso_demo.elastic_net_cost(
            c_new, w_new, t, Q, theta, lamb1, lamb2, n_dofs, include_intercept)
        df = lasso_demo.elastic_net_cost_der(
            c_new, w_new, t, Q, theta, lamb2, n_dofs, include_intercept)
        return f, df

    xopt = np.vstack((c, w))
    theta, residual, _ = lasso.multi_task_weighted_elastic_net(
        X, Q, Xdot2, alpha=a, rho=r)

    for i in range(iter_max):
        theta, c, w, p = lasso_demo.prune_params(
            c, w, theta, n_dofs, include_intercept)
        xopt = np.vstack((c, w))
        bfgs_options = {'maxiter': 100*p}
        # update RBF weights
        result = opt.minimize(f_opt, xopt, jac=True,
                              method="BFGS", options=bfgs_options)
        xopt = result.x
        c = xopt[:n_dofs*p]
        w = xopt[n_dofs*p:]
        X, Xdot2 = lasso_demo.stack_regressors(
            c, w, t, n_dofs, include_intercept)

        # perform lasso
        res_last = residual
        theta, residual, _ = lasso.multi_task_weighted_elastic_net(
            X, Q, Xdot2, alpha=a, rho=r)
        # shrink the regularizers
        # to scale lasso throughout iterations
        a /= (np.linalg.norm(res_last, 'fro') /
              np.linalg.norm(residual, 'fro'))**2
        lamb1 = 2*N*a*r
        lamb2 = N*a*(1-r)

    theta, c, w, p = lasso_demo.prune_params(
        c, w, theta, n_dofs, include_intercept)
    X, Xdot2 = lasso_demo.stack_regressors(c, w, t, n_dofs, include_intercept)
    return X, theta, c, w, a, r, Xdot2


def train_multi_demo_sparse_weights(args, p=10, plot_regr=False, path=False, examples=None, save=False, verbose=False):
    ''' Here we have multiple demonstrations '''
    import multi_dof_lasso as lasso
    joint_dict, ball_dict = serve.run_serve_demo(args)
    idx_move = joint_dict['idx_move']
    # hacky, assuming they are all 1 sec long
    Q = np.zeros((500*7, len(examples)))
    intercepts = np.zeros((7, len(examples)))
    for i, ex in enumerate(examples):
        idx = np.arange(idx_move[0, ex], idx_move[1, ex])
        q = joint_dict['x'][idx, :]
        intercepts[:, i] = np.mean(q, axis=0)
        q_flat = q.T.flatten()
        Q[:, i] = q_flat
        t = joint_dict['t'][idx]  # assumed the same for each ex
    t -= t[0]

    stacked_intercepts = np.zeros((500*7, len(examples)))
    for i, ex in enumerate(examples):
        m = np.tile(intercepts[:, i][:, np.newaxis], (1, 500))
        stacked_intercepts[:, i] = m.flatten()
    Q_zero_mean = Q - stacked_intercepts
    X, theta, c, w, a, r, Xdot2 = iter_multi_demo_lasso(t, Q_zero_mean, p)

    if verbose:
        for i, ex in enumerate(examples):
            print 'Ex:', i
            print 'Res. norm:', np.linalg.norm(
                Q_zero_mean[:, i] - np.dot(X, theta[:, i]))
            print 'Acc. norm:', np.linalg.norm(
                np.dot(Xdot2, theta[:, i]))
            print 'No. params:', theta[:, i].size

    if path:
        _, _, theta_path = lasso.multi_task_weighted_elastic_net(
            X, Q_zero_mean, Xdot2, alpha=a, rho=r, path=True)
        order = find_ordering_from_path(theta_path, plot_regr)
    else:
        order = None

    if plot_regr:
        plot_multi_regression(X, theta, Q, intercepts)

    # last one params are the intercepts!
    if save:
        dump_json_multi_demo_regression_obj(
            c, w, theta, intercepts, examples=examples, order=order, args=args)

    return


def train_multi_dof_sparse_weights(args, p=10, plot_regr=False, path=False, ex=0, save=False, verbose=False):
    ''' Train for only one demonstration but multiple dofs'''

    import multi_dof_lasso as lasso
    joint_dict, ball_dict = serve.run_serve_demo(args)

    # train multi task elastic net
    idx_move = joint_dict['idx_move']
    idx = np.arange(idx_move[0, ex], idx_move[1, ex])
    q = joint_dict['x'][idx, :]
    t = joint_dict['t'][idx]
    t -= t[0]
    c = np.linspace(t[0], t[-1], p)  # centers
    w = 0.1 * np.ones((p,))  # widths
    X = basis.create_gauss_regressor(c, w, time=t)
    C, Xdot2 = basis.create_acc_weight_mat(c, w, t, include_intercept=True)  #

    # theta, res = l2_pen_regr(X, q, C)
    # theta, res = multi_task_lasso(X, q)
    # theta, res = multi_task_elastic_net(X, q)
    # theta, res = multi_task_weighted_elastic_net(X, q, Xdot2)
    X, theta, c, w, a, r, Xdot2 = iter_multi_dof_lasso(t, q, p)

    if verbose:
        print 'Res. norm:', np.linalg.norm(q - np.dot(X, theta), 'fro')
        print 'Acc. norm:', np.linalg.norm(np.dot(Xdot2, theta), 'fro')
        print 'No. params:', theta.size

    if path:
        _, _, theta_path = lasso.multi_task_weighted_elastic_net(
            X, q, Xdot2, alpha=a, rho=r, path=True)
        order = find_ordering_from_path(theta_path)
    else:
        order = None

    if plot_regr:
        plot_regression(X, theta, q)

    # last one params are the intercepts!
    if save:
        filename = "rbf_" + str(ex) + ".json"
        dump_json_regression_obj(c, w, theta, file_name=filename, order=order)


def find_ordering_from_path(theta_path, plot):
    ''' Find the ordering of nonzero theta coefficients'''
    # print theta_path.shape
    path = theta_path[:, :, 0]  # take first matrix
    order = list()
    for i in range(path.shape[0]):
        if np.count_nonzero(path[i, :]) > 0:
            idx = np.nonzero(path[i, :])[0]
            for j in idx:
                if j not in order:
                    order.append(j)

    if plot:
        idx = np.arange(0, path.shape[0], 10)
        path_decim = path[idx, :]
        xx = np.sum(np.abs(path_decim), axis=1)
        xx /= xx[-1]
        plt.plot(xx, path_decim)
        ymin, ymax = plt.ylim()
        plt.vlines(xx, ymin, ymax, linestyle='dashed')
        plt.xlabel('$|$coef$|$ / $\max|$coef$|$')
        plt.ylabel('Coefficients')
        plt.title('Elastic Net Path')
        plt.axis('tight')

    return order


if __name__ == '__main__':
    date = '16.11.18'
    args = serve.create_default_args(date)
    args.ball_file = None
    args.num_examples = 15
    args.plot = False
    args.date = date
    examples = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    time_init = time.time()
    # train_multi_dof_sparse_weights(
    #     args, plot_regr=True, ex=18, save=False, p=500, path=False, verbose=True)
    # train_l2_reg_regr(args, plot_regr=False, ex=18,
    #                  save=False, p=10, verbose=True)
    train_multi_demo_sparse_weights(
        args, p=200, plot_regr=True, examples=examples, save=True, path=False, verbose=True)
    print 'Elapsed time:', time.time() - time_init
