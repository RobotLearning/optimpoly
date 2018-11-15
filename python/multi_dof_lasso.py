'''
Multi-task Lasso and Elastic net across dofs (degrees of freedom) of robot
are stored here
'''

import numpy as np
from sklearn import linear_model as lm
import basis_fnc as basis


def multi_task_lasso(X, q, cv=False, alpha=0.002):
    '''
    Multi Task Lasso with dimensions forced to share features
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
        #l1_ratio_list = 1-np.exp(-np.arange(1, 10)/2.0)
        clf = lm.MultiTaskElasticNetCV(l1_ratio=l1_ratio_list, eps=1e-3, n_alphas=100, alphas=None,
                                       fit_intercept=False, cv=3, verbose=True, n_jobs=-1)
    else:
        clf = lm.MultiTaskElasticNet(
            alpha=alpha, l1_ratio=l1_ratio, fit_intercept=False)
    clf.fit(X, q)
    theta = clf.coef_.T
    res = q - np.dot(X, theta)
    return theta, res


def multi_task_weighted_elastic_net(X, q, L, alpha=0.001, rho=0.99999, path=False):
    '''
    Since sklearn does not support weighted version(for us weighted L2 regularization) of Elastic Net, we transform multitask elastic net to multitask lasso
    solve it and then transform the solution back to elastic net

    L is the accelerations that composes the weight matrix, i.e. W = L'* L
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

    if path:
        res = clf.path(X_bar, q_bar, l1_ratio=0.5, n_alphas=100)
        coefs = res[1]
        theta_path = coefs.T * mult
    else:
        theta_path = []

    return theta, res, theta_path


def elastic_net_cost(c, w, t, q, theta, lamb1, lamb2):

    import basis_fnc as basis
    X = basis.create_gauss_regressor(c, w, t)
    _, Xdot2 = basis.create_acc_weight_mat(c, w, t)
    res = q - np.dot(X, theta)
    cost = np.linalg.norm(res, 'fro')**2

    theta_21_norm = np.sum(np.sqrt(np.sum(theta*theta, axis=1)))
    l2_acc_pen = np.linalg.norm(np.dot(Xdot2, theta), 'fro')**2
    cost += lamb1 * theta_21_norm
    cost += lamb2 * l2_acc_pen

    return cost


def elastic_net_cost_der(c, w, t, q, theta, lamb2):
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
