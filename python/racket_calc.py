import numpy as np

# Kinematics related functions here
# Calculate the racket orientation based on quaternion


def calcRacketOrientation(cartOrient):

    # quaternion transformation of -pi/2 from endeff to racket
    pi = np.pi
    rot = np.array([np.cos(pi/4), -np.sin(pi/4), 0, 0])
    racketOrient = mult2Quat(cartOrient, rot)
    return racketOrient


def mult2Quat(q1, q2):
    '''Function to multiply two quaternions
    '''
    q3 = np.zeros((4, 1))
    q3[0] = q1[0]*q2[0] - np.dot(q1[1:4], q2[1:4])  # q1[1:4]*q2[1:4]
    q3[1] = q1[1]*q2[0] + q1[0]*q2[1] + q1[2]*q2[3] - q1[3]*q2[2]
    q3[2] = q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1]
    q3[3] = q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0]
    return q3


def quat2Rot(q):
    ''' Computes the rotation matrices using the quaternions
    '''
    R = np.zeros((3, 3, q.shape[1]))
    R[0, 0, :] = -1.0 + 2.0*q[0, :]**2 + 2.0*q[1, :]**2
    R[1, 1, :] = -1.0 + 2.0*q[0, :]**2 + 2.0*q[2, :]**2
    R[2, 2, :] = -1.0 + 2.0*q[0, :]**2 + 2.0*q[3, :]**2

    R[0, 1, :] = 2.0*(q[1, :]*q[2, :] - q[0, :]*q[3, :])
    R[0, 2, :] = 2.0*(q[1, :]*q[3, :] + q[0, :]*q[2, :])
    R[1, 0, :] = 2.0*(q[1, :]*q[2, :] + q[0, :]*q[3, :])
    R[1, 2, :] = 2.0*(q[2, :]*q[3, :] - q[0, :]*q[1, :])
    R[2, 0, :] = 2.0*(q[1, :]*q[3, :] - q[0, :]*q[2, :])
    R[2, 1, :] = 2.0*(q[2, :]*q[3, :] + q[0, :]*q[1, :])

    return R


def rot2Quat(R):
    ''' Computes the quaternian from its rotation matrix
    '''
    R = np.squeeze(R)
    T = 1.0 + R[0, 0] + R[1, 1] + R[2, 2]

    if T > 0.00000001:

        S = 0.5 / np.sqrt(T)
        qw = 0.25 / S
        qx = (R[2, 1] - R[1, 2]) * S
        qy = (R[0, 2] - R[2, 0]) * S
        qz = (R[1, 0] - R[0, 1]) * S
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        S = np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2]) * 2
        qx = 0.25 * S
        qy = (R[0, 1] + R[1, 0]) / S
        qz = (R[0, 2] + R[2, 0]) / S
        qw = (R[1, 2] - R[2, 1]) / S
    elif R[1, 1] > R[2, 2]:
        S = np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2]) * 2
        qx = (R[0, 1] + R[1, 0]) / S
        qy = 0.25 * S
        qz = (R[1, 2] + R[2, 1]) / S
        qw = (R[0, 2] - R[2, 0]) / S
    else:
        S = np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1]) * 2
        qx = (R[0, 2] + R[2, 0]) / S
        qy = (R[1, 2] + R[2, 1]) / S
        qz = 0.25 * S
        qw = (R[0, 1] - R[1, 0]) / S

    q = np.zeros((4,))
    q[0] = qw
    q[1] = qx
    q[2] = qy
    q[3] = qz
    return q
