%% Script that loads recorded kinesthetic teach in data
% Files with 
% FIRST DATASET:
% Very good: 1, 7, 10, 16, 20
% Good: 5, 6, 9, 17, 18
%
% SECOND DATASET:
% Good: 2,4,6,9,10,12,15,16,18,19,20,21

clc; clear; close all;

file = '../data/10.6.18/joints.txt';
q = dlmread(file);
N = size(q,1);
dt = 0.002;
t = dt * (1:N);
figure;
for i = 1:7
    subplot(7,1,i);
    plot(t,q(:,i));
end

%% Find the prominent fast moving segments

qd = diff(q)./dt;
vels = sqrt(sum(qd.*qd,2));
[sorted_vels,idx_vels] = sort(vels,'descend');

size_clusters = 21; % FIRST DATASET: 20; SECOND: 21
idx_clusters = [idx_vels(1)];

idx = 2;
thresh = 1000; % 2 seconds difference minimum betw. demonstr.
while length(idx_clusters) < size_clusters
   
    diff_max = min(abs(idx_vels(idx) - idx_clusters));
    if diff_max > thresh
    
        idx_clusters = [idx_vels(idx), idx_clusters];
    end
    idx = idx + 1;
end

% sort the points and find points of low velocity around them
clusters = sort(idx_clusters);
low_vel_idxs = zeros(2,size_clusters);
low_vel_thresh = 1e-3;
max_duration_movement = 1.0; % seconds
idx_max_move = max_duration_movement/dt;
for i = 1:size_clusters
    % find the first index below cluster high vel idx
    % where the vel drops below thresh
    idx_low_pt = clusters(i);
    iter = 1;
    while (vels(idx_low_pt) > low_vel_thresh && iter < idx_max_move/2)
        idx_low_pt = idx_low_pt - 1;
        iter = iter + 1;
    end
    low_vel_idxs(1,i) = idx_low_pt;
    
    % find the first index above cluster idx
    idx_high_pt = clusters(i);
    iter = 1;
    while (vels(idx_high_pt) > low_vel_thresh && iter < idx_max_move/2)
        idx_high_pt = idx_high_pt + 1;
        iter = iter + 1;
    end
    low_vel_idxs(2,i) = idx_high_pt;
end

%% Plot the motions

examples = 10; %[2,4,6,9,10]; % FIRST: [1;7;10;5;9];
num_examples = length(examples);
q_train = cell(1,num_examples);
t_train = cell(1,num_examples);
for j = 1:length(examples)
    idx_plot = low_vel_idxs(1,examples(j)):low_vel_idxs(2,examples(j));
    q_plot = q(idx_plot,:);
    q_train{j} = q_plot;
    t_plot = t(idx_plot)';
    t_train{j} = t_plot;
    %%{
    figure('Name',['Movement #', num2str(examples(j))]);
    for i = 1:7
        subplot(7,1,i);
        plot(t_plot,q_plot(:,i));
    end
    %}
end
%{
figure('Name','Movement 1');
for i = 1:7
    subplot(7,1,i);
    plot(t_train{1},q_train{1}(:,i));
end
%}
%% Train DMPs and evolve one

n_bf = 10;
tau = 1.0; % speed up/slow down movement
dmps = trainMultiDMPs(t_train,q_train,'d',n_bf);
dmps(1).can.tau = tau;
t_evolve = 1.0/tau;
N_evolve = t_evolve/dt;
q_dmp = zeros(7,N_evolve);
ref = zeros(2*7,N_evolve);
for i = 1:7
    dmps(i).y0 = [q_train{1}(1,i);0;0];
    [x,y] = dmps(i).evolve(N_evolve);
    dmps(i).resetStates();
    q_dmp(i,:) = y(1,:);
    ref(i,:) = y(1,:)';
    ref(i+7,:) = y(2,:)';
end
t_dmp = dt*(1:N_evolve);
figure('Name','DMP evolve');
for i = 1:7
    subplot(7,1,i);
    plot(t_dmp,q_dmp(i,:));
end

%% Visualize DMP in Cartesian space

% draw Cartesian output
wam = init_wam();
y_des_cart = wam.kinematics(ref);
figure;
hold on;
grid on;
plot3(y_des_cart(1,:),y_des_cart(2,:),y_des_cart(3,:),'r--');
xlabel('x'); ylabel('y'); zlabel('z');
drawTimeIter = 40;
tLabel = t(1:drawTimeIter:N_evolve);
precision = 4;
tLabelCell = num2cell(tLabel,precision);
for i = 1:length(tLabelCell)
    tLabelCell{i} = num2str(tLabelCell{i});
end
% annotate some of the ball positions
xDraw = y_des_cart(1,1:drawTimeIter:end);
yDraw = y_des_cart(2,1:drawTimeIter:end);
zDraw = y_des_cart(3,1:drawTimeIter:end);
text(xDraw,yDraw,zDraw,tLabelCell)
scatter3(xDraw,yDraw,zDraw,20,'b','*');
hold off;

%% Save DMPs in JSON format
%{
dmp_save = struct();
dmp_save.ndofs = 7;
dmp_save.alpha = 25.0;
dmp_save.beta = 6.25;
for i = 1:7
    dmp_save.joints(i).ID = i;
    dmp_save.joints(i).init_pos = dmps(i).y0(1);
    dmp_save.joints(i).weights = dmps(i).w;
    dmp_save.joints(i).heights = dmps(1).can.h;
    dmp_save.joints(i).centers = dmps(1).can.c;
    dmp_save.joints(i).goal = dmps(i).goal;
end
json_txt = jsonencode(dmp_save);
fid = fopen('../json/dmp10.json','wt');
fwrite(fid, json_txt, 'char');
fclose(fid);
%}