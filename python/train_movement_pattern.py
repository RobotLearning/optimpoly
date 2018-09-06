'''
Train movement pattern from examples

'''
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
    plt.show()


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


def train_multi_sparse_weights(args, plot_regr=False, examples=None, save=False):
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
    X, theta, c, w = iter_multi_demo_lasso(t, Q)
    return


def train_sparse_weights(args, plot_regr=False, ex=0, p=10, save=False):
    ''' Train for only one demonstration '''
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

    theta, c, w, _ = prune_params(theta, c, w)
    # last one params are the intercepts!
    if save:
        filename = "rbf_" + str(ex) + ".json"
        dump_json_regression_obj(c, w, theta, file_name=filename)


def lasso_cost(c, w, t, q, theta, lamb1, lamb2):

    import basis_fnc as basis
    X = basis.create_gauss_regressor(c, w, t)
    C, Xdot2 = basis.create_acc_weight_mat(c, w, t)
    res = q - np.dot(X, theta)
    cost = np.linalg.norm(res, 'fro')**2

    theta_21_norm = np.sum(np.sqrt(np.sum(theta*theta, axis=1)))
    l2_acc_pen = np.linalg.norm(np.dot(Xdot2, theta), 'fro')**2
    cost += lamb1 * theta_21_norm
    cost += lamb2 * l2_acc_pen

    return cost


def lasso_cost_der(c, w, t, q, theta, lamb2):
    import basis_fnc as basis
    X = basis.create_gauss_regressor(c, w, t)
    _, Xdot2 = basis.create_acc_weight_mat(c, w, t)
    res = q - np.dot(X, theta)
    M = basis.create_gauss_regressor_der(
        c, w, t, include_intercept=True, der='c')
    Mdot2 = basis.create_acc_weight_mat_der(
        c, w, t, include_intercept=True, der='c')
    grad_c = -2 * np.diag(np.dot(np.dot(theta, res.T), M))
    grad_c += 2*lamb2 * \
        np.sum(np.dot(np.dot(Xdot2, theta), theta.T)*Mdot2, axis=0)
    M = basis.create_gauss_regressor_der(
        c, w, t, include_intercept=True, der='w')
    Mdot2 = basis.create_acc_weight_mat_der(
        c, w, t, include_intercept=True, der='w')
    grad_w = -2 * np.diag(np.dot(np.dot(theta, res.T), M))
    grad_w += 2*lamb2 * \
        np.sum(np.dot(np.dot(Xdot2, theta), theta.T)*Mdot2, axis=0)
    return np.hstack((grad_c[:-1], grad_w[:-1]))


def prune_params(theta, c, w):  # acts globally in the function
    idx_non = np.nonzero(theta[:, 0])[0]  # only nonzero entries
    p_new = len(idx_non)-1
    print 'Lasso kept', p_new, 'solutions!'  # intercept is static
    theta_new = theta[idx_non, :]
    c_new = c[idx_non[:-1]]
    w_new = w[idx_non[:-1]]
    return theta_new, c_new, w_new, p_new


def iter_multi_dof_lasso(t, q, p):
    '''
    Multi-Task grouping is performed across degrees of freedom (joints).
    Iterative MultiTaskElasticNet with nonlinear optimization (BFGS) to update BOTH the RBF
    parameters as well as the regression parameters.
    '''
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
        f = lasso_cost(c_opt, w_opt, t, q, theta, lamb1, lamb2)
        df = lasso_cost_der(c_opt, w_opt, t, q, theta, lamb2)
        return f, df

    xopt = np.hstack((c, w))
    theta, residual = multi_task_weighted_elastic_net(
        X, q, Xdot2, alpha=a, rho=r)

    for i in range(iter_max):
        theta, c, w, p = prune_params(theta, c, w)
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
        theta, residual = multi_task_weighted_elastic_net(
            X, q, Xdot2, alpha=a, rho=r)
        # shrink the regularizers
        # to scale lasso throughout iterations
        a /= (np.linalg.norm(res_last, 'fro') /
              np.linalg.norm(residual, 'fro'))**2
        lamb1 = 2*N*a*r
        lamb2 = N*a*(1-r)

    return X, theta, c, w


def iter_multi_demo_lasso(t, Q):
    '''
    Multi-Task grouping is performed across multiple demonstrations!
    Iterative MultiTaskElasticNet with nonlinear optimization (BFGS) to update BOTH the RBF
    parameters as well as the regression parameters.
    In the case of Barrett WAM there are 7x more RBF parameters than iter_multi_dof_lasso!

    Q is a IxJxK 3d-array
    T is a IxK matrix
    where I = # of time points
          J = # of dofs
          K = # of demos
    '''
    # initialize the iteration
    p = 50  # number of params. max
    n_dofs = 7  # degrees of freedom
    n_tp = Q.shape[0]/n_dofs  # number of time points
    n_demos = Q.shape[1]
    X = np.zeros((n_dofs * n_tp, p+1))
    Xdot2 = np.zeros((n_dofs * n_tp, p+1))

    # initialize big X
    c = np.linspace(t[0], t[-1], p) + 0.01 * np.random.randn(p)
    w = 0.1 * np.ones((p,)) + 0.01 * np.random.rand(p)
    for j in range(n_dofs):
        v = j*n_tp + np.arange(0, n_tp, 1)
        X[v, :] = basis.create_gauss_regressor(c, w, t)
        _, M = basis.create_acc_weight_mat(c, w, t)
        Xdot2[v, :] = M
    c = np.tile(c, (n_dofs,))
    w = np.tile(w, (n_dofs,))
    iter_max = 3
    a = 0.001  # alpha
    r = 0.99999  # rho
    N = n_tp * n_dofs
    lamb1 = 2*N*a*r
    lamb2 = N*a*(1-r)

    def f_opt(x):
        X_new = np.zeros((n_dofs * n_tp, p+1))
        C = np.zeros((p+1, p+1))
        for j in range(n_dofs):
            v_opt = j*n_tp + np.arange(0, n_tp, 1)
            c_opt = x[j*p:(j+1)*p]
            w_opt = x[(n_dofs+j)*p:(n_dofs+j+1)*p]
            X_new[v_opt, :] = basis.create_gauss_regressor(c_opt, w_opt, t)
            C += basis.create_acc_weight_mat(c_opt, w_opt, t)[0]
        res = Q - np.dot(X_new, theta)
        cost = np.linalg.norm(res, 'fro')**2
        theta_21_norm = np.sum(np.linalg.norm(theta, axis=1))
        l2_acc_pen = np.linalg.norm(
            np.dot(theta.T, np.dot(C, theta)), 'fro')**2
        cost += lamb1 * theta_21_norm
        cost += lamb2 * l2_acc_pen
        return cost

    def prune_params(params):  # acts globally in the function
        idx_non = np.nonzero(params[:, 0])[0]  # only nonzero entries
        p_new = len(idx_non)-1
        print 'Lasso kept', p_new, 'solutions!'  # intercept is static
        theta_new = params[idx_non, :]
        c_new = np.zeros((p_new*n_dofs,))
        w_new = np.zeros((p_new*n_dofs,))
        for j in range(n_dofs):
            v_opt = c[j*p:(j+1)*p]
            c_new[j*p_new:(j+1)*p_new] = v_opt[idx_non[:-1]]
            v_opt = w[j*p:(j+1)*p]
            w_new[j*p_new:(j+1)*p_new] = v_opt[idx_non[:-1]]
        return theta_new, c_new, w_new, p_new

    xopt = np.vstack((c, w))
    theta, residual = multi_task_weighted_elastic_net(
        X, Q, Xdot2, alpha=a, rho=r)

    for i in range(iter_max):
        theta, c, w, p = prune_params(theta)
        # update RBF weights
        result = opt.minimize(f_opt, xopt, method="BFGS")
        xopt = result.x

        for j in range(n_dofs):
            v = j*n_tp + np.arange(0, n_tp, 1)
            cj = xopt[j*p:(j+1)*p]
            wj = xopt[(n_dofs+j)*p:(n_dofs+j+1)*p]
            X[v, :] = basis.create_gauss_regressor(cj, wj, t)
            _, M = basis.create_acc_weight_mat(cj, wj, t)
            Xdot2[v, :] = M

        # perform lasso
        res_last = residual
        theta, residual = multi_task_weighted_elastic_net(
            X, Q, Xdot2, alpha=a, rho=r)
        # shrink the regularizers
        # to scale lasso throughout iterations
        a /= (np.linalg.norm(res_last)/np.linalg.norm(residual))**2
        lamb1 = 2*N*a*r
        lamb2 = N*a*(1-r)

    return X, theta, c, w


if __name__ == '__main__':
    date = '10.6.18'
    args = serve.create_default_args(date)
    args.ball_file = None
    args.num_examples = 22
    args.plot = False
    args.date = '10.6.18'
    #examples = [14, 15, 16, 17, 18]
    time_init = time.time()
    train_sparse_weights(args, plot_regr=True, ex=14, save=True, p=50)
    print 'Elapsed time:', time.time() - time_init
    #train_multi_sparse_weights(args, False, examples=examples, save=False)
