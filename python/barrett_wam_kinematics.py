# Barrett WAM forward kinematics
# used to show trajectories in cartesian space
#
# Function taken from SL:
# shared/barrett/math/LInfo_declare.h
# shared/barrett/math/LInfo_math.h
#
# these are called from the kinematics method of Barrett WAM class
#
# \param[out]    Xaxis   : array of rotation axes (z)
# \param[out]    Xorigin : array of coord.sys. origin vectors
# \param[out]    Xlink   : array of link position
# \param[out]    Amats   : homogeneous transformation matrices of each link
#

import numpy as np
import collections as col


def power(a, n):
    return a**n


def barrett_wam_kinematics(q):

    NDOF = 7

    # definitions
    ZSFE = 0.346  # !< z height of SAA axis above ground
    ZHR = 0.505  # !< length of upper arm until 4.5cm before elbow link
    YEB = 0.045  # !< elbow y offset
    ZEB = 0.045  # !< elbow z offset
    YWR = -0.045  # !< elbow y offset (back to forewarm)
    ZWR = 0.045  # !< elbow z offset (back to forearm)
    ZWFE = 0.255  # !< forearm length (minus 4.5cm)

    # extract parameters
    eff = col.namedtuple('struct', 'm mcm x a')
    basec = col.namedtuple('struct', 'x xd xdd')
    baseo = col.namedtuple('struct', 'q qd qdd ad add')
    eff.m = 0.0
    eff.mcm = np.zeros(3)
    eff.mcm[0] = 0.0
    eff.mcm[1] = 0.0
    eff.mcm[2] = 0.0
    eff.x = np.zeros(3)
    eff.x[0] = 0.0
    eff.x[1] = 0.0
    eff.x[2] = 0.30
    eff.a = np.zeros(3)
    eff.a[0] = 0.0
    eff.a[1] = 0.0
    eff.a[2] = 0.0
    basec.x = np.zeros((3, 1))
    basec.xd = np.zeros((3, 1))
    basec.xdd = np.zeros((3, 1))
    baseo.q = np.array([0, 1.0, 0, 0])
    baseo.qd = np.zeros((4, 1))
    baseo.qdd = np.zeros((4, 1))
    baseo.ad = np.zeros((3, 1))
    baseo.add = np.zeros((3, 1))

    # np.sine and conp.sine precomputation
    ss1th = np.sin(q[0])
    cs1th = np.cos(q[0])
    ss2th = np.sin(q[1])
    cs2th = np.cos(q[1])
    ss3th = np.sin(q[2])
    cs3th = np.cos(q[2])
    ss4th = np.sin(q[3])
    cs4th = np.cos(q[3])
    ss5th = np.sin(q[4])
    cs5th = np.cos(q[4])
    ss6th = np.sin(q[5])
    cs6th = np.cos(q[5])
    ss7th = np.sin(q[6])
    cs7th = np.cos(q[6])

    # endeffector orientations

    rseff1a1 = np.sin(eff.a[0])
    rceff1a1 = np.cos(eff.a[0])
    rseff1a2 = np.sin(eff.a[1])
    rceff1a2 = np.cos(eff.a[1])
    rseff1a3 = np.sin(eff.a[2])
    rceff1a3 = np.cos(eff.a[2])

    # Calculations are done here

    Hi00 = np.zeros((4, 4))
    Hi00[3][3] = 1
    Hi01 = np.zeros((4, 4))
    Hi01[3][3] = 1
    Hi01[2][2] = 1
    Hi12 = np.zeros((4, 4))
    Hi12[3][3] = 1
    Hi12[0][2] = -1
    Hi23 = np.zeros((4, 4))
    Hi23[3][3] = 1
    Hi23[0][2] = 1
    Hi34 = np.zeros((4, 4))
    Hi34[3][3] = 1
    Hi34[0][2] = -1
    Hi45 = np.zeros((4, 4))
    Hi45[3][3] = 1
    Hi45[0][2] = 1
    Hi56 = np.zeros((4, 4))
    Hi56[3][3] = 1
    Hi56[0][2] = -1
    Hi67 = np.zeros((4, 4))
    Hi67[3][3] = 1
    Hi67[0][2] = 1
    Hi78 = np.zeros((4, 4))
    Hi78[3][3] = 1

    # inverse homogeneous rotation matrices
    Hi00[0][0] = -1 + 2*power(baseo.q[0], 2) + 2*power(baseo.q[1], 2)
    Hi00[0][1] = 2*(baseo.q[1]*baseo.q[2] - baseo.q[0]*baseo.q[3])
    Hi00[0][2] = 2*(baseo.q[0]*baseo.q[2] + baseo.q[1]*baseo.q[3])
    Hi00[0][3] = basec.x[0]
    Hi00[1][0] = 2*(baseo.q[1]*baseo.q[2] + baseo.q[0]*baseo.q[3])
    Hi00[1][1] = -1 + 2*power(baseo.q[0], 2) + 2*power(baseo.q[2], 2)
    Hi00[1][2] = 2*(-(baseo.q[0]*baseo.q[1]) + baseo.q[2]*baseo.q[3])
    Hi00[1][3] = basec.x[1]
    Hi00[2][0] = 2*(-(baseo.q[0]*baseo.q[2]) + baseo.q[1]*baseo.q[3])
    Hi00[2][1] = 2*(baseo.q[0]*baseo.q[1] + baseo.q[2]*baseo.q[3])
    Hi00[2][2] = -1 + 2*power(baseo.q[0], 2) + 2*power(baseo.q[3], 2)
    Hi00[2][3] = basec.x[2]
    Hi01[0][0] = cs1th
    Hi01[0][1] = -ss1th
    Hi01[1][0] = ss1th
    Hi01[1][1] = cs1th
    Hi01[2][3] = ZSFE
    Hi12[1][0] = ss2th
    Hi12[1][1] = cs2th
    Hi12[2][0] = cs2th
    Hi12[2][1] = -ss2th
    Hi23[0][3] = ZHR
    Hi23[1][0] = ss3th
    Hi23[1][1] = cs3th
    Hi23[2][0] = -cs3th
    Hi23[2][1] = ss3th
    Hi34[1][0] = ss4th
    Hi34[1][1] = cs4th
    Hi34[1][3] = YEB
    Hi34[2][0] = cs4th
    Hi34[2][1] = -ss4th
    Hi34[2][3] = ZEB
    Hi45[0][3] = ZWR
    Hi45[1][0] = ss5th
    Hi45[1][1] = cs5th
    Hi45[1][3] = YWR
    Hi45[2][0] = -cs5th
    Hi45[2][1] = ss5th
    Hi56[1][0] = ss6th
    Hi56[1][1] = cs6th
    Hi56[2][0] = cs6th
    Hi56[2][1] = -ss6th
    Hi56[2][3] = ZWFE
    Hi67[1][0] = ss7th
    Hi67[1][1] = cs7th
    Hi67[2][0] = -cs7th
    Hi67[2][1] = ss7th
    Hi78[0][0] = rceff1a2*rceff1a3
    Hi78[0][1] = -(rceff1a2*rseff1a3)
    Hi78[0][2] = rseff1a2
    Hi78[0][3] = eff.x[0]
    Hi78[1][0] = rceff1a3*rseff1a1*rseff1a2 + rceff1a1*rseff1a3
    Hi78[1][1] = rceff1a1*rceff1a3 - rseff1a1*rseff1a2*rseff1a3
    Hi78[1][2] = -(rceff1a2*rseff1a1)
    Hi78[1][3] = eff.x[1]
    Hi78[2][0] = -(rceff1a1*rceff1a3*rseff1a2) + rseff1a1*rseff1a3
    Hi78[2][1] = rceff1a3*rseff1a1 + rceff1a1*rseff1a2*rseff1a3
    Hi78[2][2] = rceff1a1*rceff1a2
    Hi78[2][3] = eff.x[2]

    Ai01 = np.dot(Hi00, Hi01)
    Ai02 = np.dot(Ai01, Hi12)
    Ai03 = np.dot(Ai02, Hi23)
    Ai04 = np.dot(Ai03, Hi34)
    Ai05 = np.dot(Ai04, Hi45)
    Ai06 = np.dot(Ai05, Hi56)
    Ai07 = np.dot(Ai06, Hi67)
    Ai08 = np.dot(Ai07, Hi78)

    Ahmat = np.zeros((6, 4, 4))
    Ahmat[0, :, :] = Ai02
    Ahmat[1, :, :] = Ai03
    Ahmat[2, :, :] = Ai04
    Ahmat[3, :, :] = Ai05
    Ahmat[4, :, :] = Ai07
    Ahmat[5, :, :] = Ai08

    return Ahmat
