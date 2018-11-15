'''
Testing Least Angle Regression (LARS)
'''
import sys
sys.path.append('python/')
import numpy as np
from sklearn import linear_model as lm
from sklearn import datasets
import matplotlib.pyplot as plt


def lars_lasso_demo():
    diabetes = datasets.load_diabetes()
    X = diabetes.data
    y = diabetes.target
    #print X.shape
    # >print y.shape
    print("Computing regularization path using the LARS ...")
    alphas, _, coefs = lm.lars_path(X, y, method='lasso', verbose=True)

    xx = np.sum(np.abs(coefs.T), axis=1)
    xx /= xx[-1]

    plt.plot(xx, coefs.T)
    ymin, ymax = plt.ylim()
    plt.vlines(xx, ymin, ymax, linestyle='dashed')
    plt.xlabel('|coef| / max|coef|')
    plt.ylabel('Coefficients')
    plt.title('LASSO Path')
    plt.axis('tight')
    plt.show()

    return


def lars_on_lasso():
    ''' Lars path should coincide with lasso coordinate descent solution
    on the way
    '''
    import basis_fnc as basis
    N = 50  # num of samples
    p = 10  # dim of theta
    t = np.linspace(0, 1, N)
    c = np.linspace(t[0], t[-1], p) + 0.01 * np.random.randn(p)
    w = 0.1 * np.ones((p,)) + 0.01 * np.random.rand(p)
    X = basis.create_gauss_regressor(c, w, t, include_intercept=False)
    # _, Xdot2 = basis.create_acc_weight_mat(p, t)
    p_hidden = 4  # actual params
    np.random.seed(seed=1)
    beta = np.vstack((np.random.randn(p_hidden, 1), np.zeros((p-p_hidden, 1))))
    y = np.dot(X, beta) + 0.01 * np.random.randn(N, 1)
    y = y.squeeze()

    alphas, _, coefs = lm.lars_path(X, y, method='lasso', verbose=True)
    xx = np.sum(np.abs(coefs.T), axis=1)
    xx /= xx[-1]
    print xx
    print coefs.T

    plt.plot(xx, coefs.T)
    ymin, ymax = plt.ylim()
    plt.vlines(xx, ymin, ymax, linestyle='dashed')
    plt.xlabel('|coef| / max|coef|')
    plt.ylabel('Coefficients')
    plt.title('LASSO Path')
    plt.axis('tight')
    plt.show()

    return


if __name__ == "__main__":
    # lars_lasso_demo()
    lars_on_lasso()
