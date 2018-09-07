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


def plot_multi_regression(X, theta, Q):
    ''' X is (7N x M) matrix, theta is (MxD), Q is (7N x D) '''
    Q_smooth = np.dot(X, theta)

    N = X.shape[0]/7
    D = theta.shape[1]
    t = 0.002 * np.arange(1, N+1)
    for d in range(D):
        q = Q[:, d]
        q = np.reshape(q, ((7, N)))
        q_smooth = Q_smooth[:, d]
        q_smooth = np.reshape(q_smooth, ((7, N)))
        f, axs = plt.subplots(7, 1, sharex=True)
        for i in range(7):
            axs[i].plot(t, q[i, :])
            axs[i].plot(t, q_smooth[i, :])
    plt.show(block=False)


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


def train_multi_demo_sparse_weights(args, p=10, plot_regr=False, examples=None, save=False):
    ''' Here we have multiple demonstrations '''
    joint_dict, ball_dict = serve.run_serve_demo(args)
    idx_move = joint_dict['idx_move']
    # hacky, assuming they are all 1 sec long
    Q = np.zeros((500*7, len(examples)))

    for i, ex in enumerate(examples):
        idx = np.arange(idx_move[0, ex], idx_move[1, ex])
        q = joint_dict['x'][idx, :]
        Q[:, i] = q.T.flatten()
        t = joint_dict['t'][idx]  # assumed the same for each ex
    t -= t[0]
    X, theta, c, w = iter_multi_demo_lasso(t, Q, p)
    if plot_regr:
        plot_multi_regression(X, theta, Q)

    # last one params are the intercepts!
    if save:
        for i, ex in enumerate(examples):
            filename = "rbf_" + str(ex) + ".json"
            dump_json_regression_obj(c, w, theta[:i], file_name=filename)

    return


def train_multi_dof_sparse_weights(args, plot_regr=False, ex=0, p=10, save=False):
    ''' Train for only one demonstration but multiple dofs'''
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
    X, theta, c, w = iter_multi_dof_lasso(t, q, p)
    if plot_regr:
        plot_regression(X, theta, q)

    # last one params are the intercepts!
    if save:
        filename = "rbf_" + str(ex) + ".json"
        dump_json_regression_obj(c, w, theta, file_name=filename)


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
    theta, residual = lasso.multi_task_weighted_elastic_net(
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
        theta, residual = lasso.multi_task_weighted_elastic_net(
            X, q, Xdot2, alpha=a, rho=r)
        # shrink the regularizers
        # to scale lasso throughout iterations
        a /= (np.linalg.norm(res_last, 'fro') /
              np.linalg.norm(residual, 'fro'))**2
        lamb1 = 2*N*a*r
        lamb2 = N*a*(1-r)

    theta, c, w, p = lasso.prune_params(theta, c, w)
    X = basis.create_gauss_regressor(c, w, t)
    return X, theta, c, w


def iter_multi_demo_lasso(t, Q, p):
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
    '''
    import multi_dof_lasso as lasso
    import multi_demo_lasso as lasso_demo
    # initialize the iteration
    n_dofs = 7  # degrees of freedom
    n_tp = Q.shape[0]/n_dofs  # number of time points
    n_demos = Q.shape[1]

    # initialize big X
    c = np.linspace(t[0], t[-1], p) + 0.01 * np.random.randn(p)
    w = 0.1 * np.ones((p,)) + 0.01 * np.random.rand(p)
    c = np.tile(c, (n_dofs,))
    w = np.tile(w, (n_dofs,))
    X, Xdot2 = lasso_demo.stack_regressors(c, w, t, n_dofs)
    iter_max = 3
    a = 0.001  # alpha
    r = 0.99999  # rho
    N = n_tp * n_dofs
    lamb1 = 2*N*a*r
    lamb2 = N*a*(1-r)

    def f_opt(x):
        c_new = x[:n_dofs*p]
        w_new = x[n_dofs*p:]
        f = lasso_demo.elastic_net_cost(
            c_new, w_new, t, Q, theta, lamb1, lamb2)
        df = lasso_demo.elastic_net_cost_der(
            c_new, w_new, t, Q, theta, lamb2, n_dofs)
        return f, df

    xopt = np.vstack((c, w))
    theta, residual = lasso.multi_task_weighted_elastic_net(
        X, Q, Xdot2, alpha=a, rho=r)

    bfgs_options = {'maxiter': 1000}
    for i in range(iter_max):
        theta, c, w, p = lasso_demo.prune_params(c, w, theta, n_dofs)
        xopt = np.vstack((c, w))
        # update RBF weights
        result = opt.minimize(f_opt, xopt, jac=True,
                              method="BFGS", options=bfgs_options)
        xopt = result.x
        c = xopt[:n_dofs*p]
        w = xopt[n_dofs*p:]
        X, Xdot2 = lasso_demo.stack_regressors(c, w, t, n_dofs)

        # perform lasso
        res_last = residual
        theta, residual = lasso.multi_task_weighted_elastic_net(
            X, Q, Xdot2, alpha=a, rho=r)
        # shrink the regularizers
        # to scale lasso throughout iterations
        a /= (np.linalg.norm(res_last, 'fro') /
              np.linalg.norm(residual, 'fro'))**2
        lamb1 = 2*N*a*r
        lamb2 = N*a*(1-r)

    theta, c, w, p = lasso_demo.prune_params(c, w, theta, n_dofs)
    X, _ = lasso_demo.stack_regressors(c, w, t, n_dofs)
    return X, theta, c, w


if __name__ == '__main__':
    date = '10.6.18'
    args = serve.create_default_args(date)
    args.ball_file = None
    args.num_examples = 22
    args.plot = False
    args.date = '10.6.18'
    examples = [14, 15, 16, 17, 18]
    time_init = time.time()
    # train_multi_dof_sparse_weights(args, plot_regr=True, ex=14, save=True, p=50)
    train_multi_demo_sparse_weights(
        args, p=500, plot_regr=True, examples=examples, save=False)
    print 'Elapsed time:', time.time() - time_init
