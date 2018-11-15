'''
Multi Demo version of Elastic Net/Lasso functions are stored here
where joints (dofs of robot) are stacked horizontally and in the multi-task dimension
are the different demonstrations
'''
import numpy as np
import basis_fnc as basis


def elastic_net_cost(c, w, t, q, theta, lamb1, lamb2, ndofs=7, include_intercept=True):
    ''' Unlike the lasso cost, this one contains the combined costs
    of different demonstrations and c, w are different for each dof'''
    X, Xdot2 = stack_regressors(c, w, t, ndofs, include_intercept)
    res = q - np.dot(X, theta)
    cost = np.linalg.norm(res, 'fro')**2
    # np.sum(np.linalg.norm(theta, axis=1))
    theta_21_norm = np.sum(
        np.sqrt(np.sum(theta*theta, axis=1)))
    l2_acc_pen = np.linalg.norm(np.dot(Xdot2, theta), 'fro')**2
    cost += lamb1 * theta_21_norm
    cost += lamb2 * l2_acc_pen
    return cost


def elastic_net_cost_der(c, w, t, q, theta, lamb2, n_dofs, include_intercept=True):

    X, Xdot2 = stack_regressors(c, w, t, n_dofs, include_intercept)
    res = q - np.dot(X, theta)
    M, Mdot2 = stack_der_regressors(
        c, w, t, der='c', n_dofs=n_dofs, include_intercept=include_intercept)
    if include_intercept:
        p = theta.shape[0]-1
    else:
        p = theta.shape[0]
    n_tp = len(t)
    grad_cs = np.zeros((p*n_dofs,))
    grad_ws = np.zeros((p*n_dofs,))
    for i in range(n_dofs):  # ndofs stacked horizontally
        v = i*n_tp + np.arange(0, n_tp, 1)
        grad_c = -2 * np.diag(np.dot(np.dot(theta, res[v, :].T), M[v, :]))
        grad_c += 2*lamb2 * \
            np.sum(np.dot(np.dot(Xdot2[v, :], theta),
                          theta.T)*Mdot2[v, :], axis=0)
        if include_intercept:
            grad_cs[i*p:(i+1)*p] = grad_c[:-1]
        else:
            grad_cs[i*p:(i+1)*p] = grad_c

    M, Mdot2 = stack_der_regressors(
        c, w, t, der='w', n_dofs=n_dofs, include_intercept=include_intercept)
    for i in range(n_dofs):  # ndofs stacked horizontally
        v = i*n_tp + np.arange(0, n_tp, 1)
        grad_w = -2 * np.diag(np.dot(np.dot(theta, res[v, :].T), M[v, :]))
        grad_w += 2*lamb2 * \
            np.sum(np.dot(np.dot(Xdot2[v, :], theta),
                          theta.T)*Mdot2[v, :], axis=0)
        if include_intercept:
            grad_ws[i*p:(i+1)*p] = grad_w[:-1]
        else:
            grad_ws[i*p:(i+1)*p] = grad_w
    grad = np.hstack((grad_cs, grad_ws))
    return grad


def stack_der_regressors(c, w, t, der, n_dofs, include_intercept=True):
    '''
    Form X matrix and its second time derivative Xdot2's derivatives
    w.r.t. basis function parameters c and w
    Here unlike multi_dof case, X matrices are stacked horizontally
    for each joint, and c,w, are n_dof times the normal amount
    '''
    n_tp = len(t)
    p = len(c)/n_dofs
    if include_intercept:
        X = np.zeros((n_dofs * n_tp, p+1))
        Xdot2 = np.zeros((n_dofs * n_tp, p+1))
    else:
        X = np.zeros((n_dofs * n_tp, p))
        Xdot2 = np.zeros((n_dofs * n_tp, p))

    for j in range(n_dofs):
        v = j*n_tp + np.arange(0, n_tp, 1)
        c_dof = c[j*p:(j+1)*p]
        w_dof = w[j*p:(j+1)*p]
        X[v, :] = basis.create_gauss_regressor_der(
            c_dof, w_dof, t, include_intercept, der)
        M = basis.create_acc_weight_mat_der(
            c_dof, w_dof, t, include_intercept, der)
        Xdot2[v, :] = M
    return X, Xdot2


def stack_regressors(c, w, t, n_dofs=7, include_intercept=True):
    '''
    For multi-demo grouping, form X matrix and its second derivative
    by stacking X's corresponding to different dofs
    '''
    n_tp = len(t)
    p = len(c)/n_dofs
    if include_intercept:
        X = np.zeros((n_dofs * n_tp, p+1))
        Xdot2 = np.zeros((n_dofs * n_tp, p+1))
    else:
        X = np.zeros((n_dofs * n_tp, p))
        Xdot2 = np.zeros((n_dofs * n_tp, p))

    for j in range(n_dofs):
        v = j*n_tp + np.arange(0, n_tp, 1)
        c_dof = c[j*p:(j+1)*p]
        w_dof = w[j*p:(j+1)*p]
        X[v, :] = basis.create_gauss_regressor(
            c_dof, w_dof, t, include_intercept)
        _, M = basis.create_acc_weight_mat(c_dof, w_dof, t, include_intercept)
        Xdot2[v, :] = M
    return X, Xdot2


def prune_params(c, w, params, n_dofs=7, include_intercept=True):
    idx_non = np.nonzero(params[:, 0])[0]  # only nonzero entries
    theta_new = params[idx_non, :]
    if include_intercept:
        p = params.shape[0]-1
        idx_non = idx_non[:-1]
    else:
        p = params.shape[0]
    p_new = len(idx_non)
    print 'Lasso kept', p_new, 'solutions!'  # intercept is static
    c_new = np.zeros((p_new*n_dofs,))
    w_new = np.zeros((p_new*n_dofs,))
    for j in range(n_dofs):
        v_opt = c[j*p:(j+1)*p]
        c_new[j*p_new:(j+1)*p_new] = v_opt[idx_non]
        v_opt = w[j*p:(j+1)*p]
        w_new[j*p_new:(j+1)*p_new] = v_opt[idx_non]
    return theta_new, c_new, w_new, p_new
