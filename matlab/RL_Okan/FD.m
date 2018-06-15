%% My own implementation of Finite Difference based RL
% First step is to get the same sort of learning curve
clc; clear; close all;

goal = 1.0;
reward = @(traj,goal) sum(exp(-abs(goal-traj)))/length(traj);
num_param = 10;
num_iter = 2000;
num_batch = 2*num_param;
num_steps = 100;
var_explore = 0.01;
learn_rate = 1000;
dt = 0.01;
alpha = 25;
beta = 25/4;
yin = zeros(3,1);
ax = 1.0;
tau = 1.0;
can = CAN(dt,ax,tau,num_param,1.0,'d');
dmp = DDMP(can,alpha,beta,goal,yin);
% initialize weights to zero
dmp.setWeights(zeros(1,num_param));
R = zeros(1,num_iter);
R_batch = zeros(1,num_batch);
params = zeros(num_param,num_iter+1);
X_perturb = zeros(num_batch,num_param);
y_perturb = zeros(num_batch,1);
for i = 1:num_iter
    if rem(i,100) == 0
        fprintf('Iteration %d\n',i);
    end
    
    for j = 1:num_batch
        X_perturb(j,:) = sqrt(var_explore)*randn(1,num_param);
        w = params(:,i) + X_perturb(j,:)';
        dmp.setWeights(w);
        [~,y] = dmp.evolve(num_steps);
        traj = y(1,:);
        R_batch(j) = reward(traj,goal);
    end
    y_perturb = R_batch - R(i);
    %grad = X_perturb \ y_perturb(:);
    grad = [X_perturb,ones(num_batch,1)] \ y_perturb(:);
    params(:,i+1) = params(:,i) + learn_rate * grad(1:num_param);
    
    dmp.setWeights(params(:,i+1));
    [t,y] = dmp.evolve(num_steps);
    traj = y(1,:);
    R(i) = reward(traj,goal);
end
plot(1:num_iter,R); 
%[~,y] = dmp.evolve(num_steps);
%plot(1:num_steps,y(1,:));