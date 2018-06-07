%% Generic script to initialize Barrett WAM
% alpha is the multiplier for the R penalty matrix
function [wam,PD,Q0] = init_wam()

dyn.use_cpp = true; 
dyn.nom = 1; 
dyn.act = 3; %not needed here!
N_DOFS = 7;
% Simulation Values 
% system matrices are continous
SIM.discrete = false;
% learn in cartesian space
SIM.cartesian = false;
% dimension of the x vector
SIM.dimx = 2*N_DOFS;
% dimension of the output y
SIM.dimy = 2*N_DOFS;
% dimension of the control input
SIM.dimu = N_DOFS;
% time step h 
SIM.h = 0.002; % 500 Hz recorded data
% measurement noise covariance (eps_m * eye)
SIM.eps_m = 0e-10;
% integration method
SIM.int = 'Symplectic Euler';
% reference trajectory in joint space?
SIM.jref = true;
% We observe all joints and all joint velocities
SIM.C = eye(SIM.dimy,SIM.dimx);

% form constraints
MAX_VEL = 10;
MAX_ACC = 200;
SLACK = 0.05;
CON.q.max = [2.60; 2.00; 2.80; 3.10; 1.30; 1.60; 2.20] - SLACK;
CON.q.min = [-2.60; -2.00; -2.80; -0.90; -4.80; -1.60; -2.20] + SLACK;
CON.qd.max = MAX_VEL * ones(7,1);
CON.qd.min = -MAX_VEL * ones(7,1);
CON.qdd.max = MAX_ACC * ones(7,1);
CON.qdd.min = -MAX_ACC * ones(7,1);
CON.u.max = [75; 125; 39; 30; 3; 4; 1];
CON.u.min = -CON.u.max; % ONLY Q and QD constraints are used in opt.

% initialize model
wam = BarrettWAM(CON,SIM,dyn);

% PD control defined here
PD = zeros(N_DOFS,2*N_DOFS);
PD(1,1) = -200;
PD(1,N_DOFS+1) = -7.0;
PD(2,2) = -300;
PD(2,N_DOFS+2) = -15.0;
PD(3,3) = -100;
PD(3,N_DOFS+3) = -5.0;
PD(4,4) = -50;
PD(4,N_DOFS+4) = -2.5;
PD(5,5) = -10;
PD(5,N_DOFS+5) = -0.3;
PD(6,6) = -10;
PD(6,N_DOFS+6) = -0.3;
PD(7,7) = -2.5;
PD(7,N_DOFS+7) = -0.075;

% initialize the arm with zero velocity on the right hand side
q_right = [1.0; -0.2; -0.1; 1.8; -1.57; 0.1; 0.3];
q_left = [-1.0; 0.0; 0.0; 1.5; -1.57; 0.1; 0.3];
q_centre = [0.0; 0.0; 0.0; 1.5; -1.75; 0.0; 0.0];
q0 = q_right;
qd0 = zeros(7,1);
Q0 = [q0;qd0];