'''
test transformation of weighted regression to normal regression
'''
import sys
sys.path.append('python/')
import numpy as np


def test_weighted_regr():
    '''
    Weighted regression can be solved using cholesky tranform
    and standard regression
    '''
    from numpy import eye
    from numpy import dot
    from numpy.linalg import solve
    from numpy.linalg import cholesky
    from numpy.linalg import inv
    from numpy.random import randn
    from numpy.random import rand

    n = 10
    d = 5
    X = randn(n, d)
    y = randn(n, 1)
    S = 2*eye(n) + 0.2*rand(n, n)  # inverse covar
    S = (S + S.T)/2.0
    W = 5*eye(d) + 0.2*rand(d, d)  # inverse covar of prior
    W = (W + W.T)/2.0
    beta1 = solve(dot(X.T, dot(S, X)) + W, dot(X.T, dot(S, y)))

    # compare with cholesky based soln
    M = inv(cholesky(W).T)
    Xbar = dot(cholesky(S).T, dot(X, M))
    ybar = dot(cholesky(S).T, y)
    beta2_bar = solve(dot(Xbar.T, Xbar) + eye(d), dot(Xbar.T, ybar))
    beta2 = dot(M, beta2_bar)
    assert np.allclose(beta1, beta2)


def test_basis_fnc_second_time_der():
    ''' Second derivative must match num. derivative'''
    import basis_fnc as basis
    N = 10
    c = np.linspace(1, 20, N)
    w = 10.0 * np.ones((N,))
    t = np.arange(1, 11)
    h = 1e-6
    num_der2 = np.zeros((N,))
    der2 = np.zeros((N,))
    for i in range(N):
        val_plus = basis.basis_fnc_gauss(t[i]+2*h, c[i], w[i])
        val_minus = basis.basis_fnc_gauss(t[i]-2*h, c[i], w[i])
        val = basis.basis_fnc_gauss(t[i], c[i], w[i])
        num_der2[i] = (val_plus - 2*val + val_minus) / (4*h*h)
        der2[i] = basis.basis_fnc_gauss_der2(t[i], c[i], w[i])
    assert np.allclose(num_der2, der2, atol=1e-3)


def test_elastic_net_multi_dof_gradient():
    '''
    Test the gradient of the multi_dof (only one demo!) version of the 
    Elastic Net cost with respect to
    centers and widths (i.e. parameters) by comparing to num. derivative
    '''
    import multi_dof_lasso as lasso
    N = 10
    p = 4
    d = 5  # number of columns of q
    t = np.linspace(0, 1, N)
    c = np.linspace(t[0], t[-1], p) + 0.01*np.random.randn(p)
    w = 0.1 * np.ones((p,)) + 0.01 * np.random.rand(p)
    h = 1e-4
    num_der = np.zeros((2*p,))
    theta = np.random.randn(p+1, d)
    q = np.random.randn(N, d)
    lamb1 = 1e-3
    lamb2 = 1e-3
    for i in range(p):
        c[i] += h
        f_plus = lasso.elastic_net_cost(c, w, t, q, theta, lamb1, lamb2)
        c[i] -= 2*h
        f_minus = lasso.elastic_net_cost(
            c, w, t, q, theta, lamb1, lamb2)
        num_der[i] = (f_plus - f_minus)/(2*h)
        c[i] += h
        w[i] += h
        f_plus = lasso.elastic_net_cost(
            c, w, t, q, theta, lamb1, lamb2)
        w[i] -= 2*h
        f_minus = lasso.elastic_net_cost(
            c, w, t, q, theta, lamb1, lamb2)
        num_der[p+i] = (f_plus - f_minus)/(2*h)
        w[i] += h
    exact_der = lasso.elastic_net_cost_der(c, w, t, q, theta, lamb2)
    assert np.allclose(num_der, exact_der, atol=1e-3)


def test_elastic_net_multi_demo_gradient():
    '''
    Test the gradient of the multi_demo (stacked dofs horizontally!) version of the 
    Elastic Net cost with respect to
    centers and widths (i.e. parameters) of different joints by comparing to num. derivative
    '''
    import multi_demo_lasso as lasso
    n_tps = 6  # number of timxe points
    p = 5  # number of parameters (except for intercept)
    n_dofs = 3
    n_demos = 4  # number of demonstrations
    t = np.linspace(0, 1, n_tps)
    c = np.linspace(t[0], t[-1], p) + 0.01*np.random.randn(p)
    w = 0.1 * np.ones((p,)) + 0.01 * np.random.rand(p)
    c = np.tile(c, (n_dofs,))
    w = np.tile(w, (n_dofs,))
    h = 1e-4
    num_der = np.zeros((2*p*n_dofs,))
    theta = np.random.randn(p+1, n_demos)
    q = np.random.randn(n_tps*n_dofs, n_demos)
    lamb1 = 1e-3
    lamb2 = 1e-3
    for i in range(p*n_dofs):
        c[i] += h
        f_plus = lasso.elastic_net_cost(
            c, w, t, q, theta, lamb1, lamb2, n_dofs)
        c[i] -= 2*h
        f_minus = lasso.elastic_net_cost(
            c, w, t, q, theta, lamb1, lamb2, n_dofs)
        num_der[i] = (f_plus - f_minus)/(2*h)
        c[i] += h
        w[i] += h
        f_plus = lasso.elastic_net_cost(
            c, w, t, q, theta, lamb1, lamb2, n_dofs)
        w[i] -= 2*h
        f_minus = lasso.elastic_net_cost(
            c, w, t, q, theta, lamb1, lamb2, n_dofs)
        num_der[p*n_dofs+i] = (f_plus - f_minus)/(2*h)
        w[i] += h
    exact_der = lasso.elastic_net_cost_der(c, w, t, q, theta, lamb2, n_dofs)
    assert np.allclose(num_der, exact_der, atol=1e-3)


def test_elastic_net_to_lasso_transform():
    ''' The solutions to elastic net and lasso should be identical
    after transformation
    '''
    import sklearn.linear_model as lm
    import basis_fnc as basis
    N = 50  # num of samples
    p = 10  # dim of theta
    t = np.linspace(0, 1, N)
    c = np.linspace(t[0], t[-1], p) + 0.01 * np.random.randn(p)
    w = 0.1 * np.ones((p,)) + 0.01 * np.random.rand(p)
    X = basis.create_gauss_regressor(c, w, t, include_intercept=False)
    # _, Xdot2 = basis.create_acc_weight_mat(p, t)
    p_hidden = 4  # actual params
    np.random.seed(seed=1)  # 10 passes
    beta = np.vstack((np.random.randn(p_hidden, 1), np.zeros((p-p_hidden, 1))))
    y = np.dot(X, beta) + 0.01 * np.random.randn(N, 1)
    alpha_elastic = 0.3
    ratio = 0.5

    clf = lm.ElasticNet(alpha=alpha_elastic,
                        l1_ratio=ratio, fit_intercept=False)
    clf.fit(X, y)
    beta_hat_1 = clf.coef_

    lamb1 = 2*N*alpha_elastic*ratio
    lamb2 = N*alpha_elastic*(1-ratio)
    y_bar = np.vstack((y, np.zeros((p, 1))))  # if unweighted
    # y_bar = np.vstack((y, np.zeros((N, 1))))
    mult = np.sqrt(1.0/(1+lamb2))
    X_bar = mult * np.vstack((X, np.sqrt(lamb2)*np.eye(p)))  # if unweighted
    # X_bar = mult * np.vstack((X, np.sqrt(lamb2)*Xdot2))
    lamb_bar = lamb1 * mult
    alpha_lasso = lamb_bar/(2*N)

    clf2 = lm.Lasso(alpha=alpha_lasso, fit_intercept=False)
    clf2.fit(X_bar, y_bar)
    beta_hat_2 = clf2.coef_ * mult  # transform back
    print 'Actual param:', beta.T
    print 'Elastic net est:', beta_hat_1
    print 'Lasso est:', beta_hat_2

    # print mult
    # return X
    assert np.allclose(beta_hat_1, beta_hat_2, atol=1e-3)
