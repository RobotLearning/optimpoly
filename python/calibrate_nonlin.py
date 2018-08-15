import numpy as np
import scipy.optimize as opt
import time


def est_calib_params_nonlin(pts3d,
                            pts2d,
                            distortion_dict=None,
                            intrinsic_dict=None,
                            extrinsic_dict=None,
                            num_iter_max=100,
                            debug=False):
    ''' Nonlinear estimation of calibration parameters with BFGS
    I would have liked to use scipy's curve_fit nonlinear least squares routine
    but it requires a scalar dependent variable, in our case we have both x and y!

    The parameters should probably be initialized with a linear estimation procedure
    (see svd-based - with a 2Nx12 matrix - solution of projection matrix estimation)

    --- Distortion ---
    k1-k6: distortion parameters (radial)
    p1-p2: tangential distortion

    --- Projection params ---
    R,t: rotation and translation [3 Euler angles and 3 translation params]
    intrinsic matrix: upper triangular 3x3 matrix with 1 at [2,2] position

    --- Pixels and 3d positions
    pts2d: 2 x N matrix
    pts3d: 3 x N matrix
    '''

    if distortion_dict is None:
        distortion_dict = {'k': np.zeros((6, 1)), 'p': np.zeros((2, 1))}

    if extrinsic_dict is None:
        extrinsic_dict = {'euler_angles': np.zeros(
            (3, 1)), 'trans_vector': np.zeros((3, 1))}

    if intrinsic_dict is None:
        intrinsic_dict = {'fx': 0, 'fy': 0,
                          'shear': 0, 'u0': 0, 'v0': 0, 'scale': 1}

    x0 = np.vstack((distortion_dict['k'], distortion_dict['p']))
    x0_ = np.array([intrinsic_dict['fx'],
                    intrinsic_dict['fy'],
                    intrinsic_dict['shear'],
                    intrinsic_dict['u0'],
                    intrinsic_dict['v0'],
                    intrinsic_dict['scale'],
                    extrinsic_dict['euler_angles'][0],
                    extrinsic_dict['euler_angles'][1],
                    extrinsic_dict['euler_angles'][2],
                    extrinsic_dict['trans_vector'][0],
                    extrinsic_dict['trans_vector'][1],
                    extrinsic_dict['trans_vector'][2]])
    x0 = np.hstack((x0.reshape((len(x0),)), x0_))
    options = {'maxiter': num_iter_max, 'disp': debug}
    if debug:
        call_back = print_info
    else:
        def call_back(theta): return

    num_param = len(x0)
    bounds = [(-np.inf, np.inf)]*num_param  # [None]*num_param
    # for i in range(8):  # for k and p distortion params, restrict to [-1,1]
    #    bounds[i] = (-1, 1)
    for i in range(8):
        bounds[i] = (0, 0)  # set distortion params to zero
    bounds = tuple(bounds)
    result = opt.minimize(method='L-BFGS-B', args=(pts2d, pts3d),  # opt.fmin_bfgs(
                          fun=residual, x0=x0, options=options,
                          bounds=bounds,
                          callback=call_back)  # jac=res_der
    xopt = result.x
    print(result.message)
    print('residual after opt = ', result.fun)
    distortion_dict = {'k': xopt[0:6],
                       'p': xopt[6:8]}
    intrinsic_dict = {'fx': xopt[8], 'fy': xopt[9],
                      'shear': xopt[10], 'u0': xopt[11], 'v0': xopt[12],
                      'scale': xopt[13]}
    extrinsic_dict = {'euler_angles': xopt[14:17], 'trans_vector': xopt[17:]}

    params = dict()
    params['distortion'] = distortion_dict
    params['intrinsic'] = intrinsic_dict
    params['extrinsic'] = extrinsic_dict

    return params


def print_info(theta):
    ''' Call back function to print info'''

    # time it
    if 'tic' not in print_info.__dict__:
        print_info.tic = time.time()
        print_info.iter = 1
    else:
        toc = time.time()
        print('Time elapsed (ms) :', 1000*(toc - print_info.tic))
        print_info.tic = time.time()
        print_info.iter = print_info.iter + 1

    ''' Print optim info'''
    k = theta[0:6]
    p = theta[6:8]
    P = form_extrinsic_mat(theta)
    A = form_intrinsic_mat(theta)

    print('==================ITER: {}=============='.format(print_info.iter))
    print('k: {}, \np: {}'.format(k, p))
    print('Extrinsic mat ([R,t]) = \n{}'.format(P))
    print('Intrinscit mat = \n{}'.format(A))


def form_rot_matrix(euler_angles):
    ''' Form the rotation matrix from 3 euler angles
    for some reson we have to also take transpose!!'''
    mat0 = np.eye(3)
    mat0[1, 1] = np.cos(euler_angles[0])
    mat0[1, 2] = np.sin(euler_angles[0])
    mat0[2, 1] = -mat0[1, 2]
    mat0[2, 2] = mat0[1, 1]

    mat1 = np.eye(3)
    mat1[0, 0] = np.cos(euler_angles[1])
    mat1[0, 2] = -np.sin(euler_angles[1])
    mat1[2, 0] = -mat1[0, 2]
    mat1[2, 2] = mat1[0, 0]

    mat2 = np.eye(3)
    mat2[0, 0] = np.cos(euler_angles[2])
    mat2[0, 1] = np.sin(euler_angles[2])
    mat2[1, 0] = -mat2[0, 1]
    mat2[1, 1] = mat2[0, 0]

    return np.dot(mat0, np.dot(mat1, mat2)).T


def form_extrinsic_mat(xopt):
    ''' form the matrix composed of rotation and translation (scaling is 1)'''
    euler_angles = xopt[14:17]
    trans = xopt[17:].reshape((3, 1))
    rot_mat = form_rot_matrix(euler_angles)
    return np.hstack((rot_mat, trans))


def form_intrinsic_mat(xopt):
    ''' forming the camera (intrinsic) matrix, and scaling'''
    mat = np.zeros((3, 3))
    mat[0, 0] = xopt[8]
    mat[0, 1] = xopt[10]  # shear
    mat[0, 2] = xopt[11]
    mat[1, 1] = xopt[9]
    mat[1, 2] = xopt[12]
    mat[2, 2] = 1.0
    mat = xopt[13] * mat  # scale
    return mat


def residual(theta, pts2d, pts3d):
    ''' Return residual of fitting to pixels given theta parameters'''
    res = 0
    N = pts3d.shape[1]
    P = form_extrinsic_mat(theta)
    A = form_intrinsic_mat(theta)
    k = theta[0:6]
    p = theta[6:8]
    # numpy broadcasting

    vec3s = np.dot(P, np.vstack((pts3d, np.ones((1, N)))))
    xps = vec3s[0, :]/vec3s[-1, :]
    yps = vec3s[1, :]/vec3s[-1, :]
    '''
    # ADD EFFECT OF DISTORTION MODEL
    r2s = xps**2 + yps**2
    r4s = r2s**2
    r6s = r2s*r4s
    facts = (1 + k[0]*r2s + k[1]*r4s + k[2]*r6s) / \
        (1 + k[3]*r2s + k[4]*r4s + k[5]*r6s)
    xpps = xps * facts + 2*p[0]*xps*yps + p[1]*(r2s + 2*(xps**2))
    ypps = yps * facts + p[0]*(r2s + 2*(yps**2)) + 2*p[1]*xps*yps
    '''
    xpps = xps
    ypps = yps
    pts2d_hat = np.dot(A, np.vstack((xpps, ypps, np.ones((1, N)))))
    difs = pts2d - pts2d_hat[0:-1, :]
    res = np.sum(difs * difs)
    return res


def res_der(theta, pts2d, pts3d):
    ''' Calculate derivative of the residual function'''
    # TODO:


def test_function_derivative():
    ''' Check the derivative of the residual
    by comparing with numerical derivative
    '''
    print('TODO!')


if __name__ == '__main__':
    test_function_derivative()
