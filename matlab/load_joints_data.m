%% Script that loads recorded kinesthetic teach in data
% Files with 
% FIRST DATASET:
% Very good: 1, 7, 10, 16, 20
% Good: 5, 6, 9, 17, 18
%
% SECOND DATASET:
% Good: 2,4,6,9,10,12,15,16,18,19,20,21

clc; clear; close all;

file = '../data/15.11.18/joints.txt';
ts = dlmread(file);
q = ts(:,2:end);
t = ts(:,1);
t = t - t(1);
t = 0.001 * t;
dt = 0.002;
% N = size(q,1);
% dt = 0.002;
% t = dt * (1:N);
% figure;
% for i = 1:7
%     subplot(7,1,i);
%     plot(t,q(:,i));
% end

%% Find the prominent fast moving segments

qd = diff(q)./dt;
vels = sqrt(sum(qd.*qd,2));
[sorted_vels,idx_vels] = sort(vels,'descend');

num_movements = 15; 
idx_clusters = [idx_vels(1)];

idx = 2;
thresh = 1000; % 2 seconds difference minimum betw. demonstr.
while length(idx_clusters) < num_movements
   
    diff_max = min(abs(idx_vels(idx) - idx_clusters));
    if diff_max > thresh
    
        idx_clusters = [idx_vels(idx), idx_clusters];
    end
    idx = idx + 1;
end

% sort the points and find points of low velocity around them
clusters = sort(idx_clusters);
low_vel_idxs = zeros(2,num_movements);
low_vel_thresh = 1e-3;
max_duration_movement = 1.0; % seconds
idx_max_move = max_duration_movement/dt;
for i = 1:num_movements
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

plot_examples = 1:15; 
num_examples = length(plot_examples);
q_train = cell(1,num_examples);
t_train = cell(1,num_examples);
for j = 1:length(plot_examples)
    idx_plot = low_vel_idxs(1,plot_examples(j)):low_vel_idxs(2,plot_examples(j));
    q_plot = q(idx_plot,:);
    q_train{j} = q_plot;
    t_plot = t(idx_plot)';
    t_train{j} = t_plot;
    %%{
    figure('Name',['Movement #', num2str(plot_examples(j))]);
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
% draw Cartesian output
%{
figure('Name','Training cartesian data');
hold on;
grid on;
plot3(train_cart(1,:),train_cart(2,:),train_cart(3,:),'r--');
xlabel('x'); ylabel('y'); zlabel('z');
drawTimeIter = 40;
tLabel = t_plot(1:drawTimeIter:max_duration_movement/dt);
precision = 4;
tLabelCell = num2cell(tLabel,precision);
for i = 1:length(tLabelCell)
    tLabelCell{i} = num2str(tLabelCell{i});
end
% annotate some of the ball positions
xDraw = train_cart(1,1:drawTimeIter:end);
yDraw = train_cart(2,1:drawTimeIter:end);
zDraw = train_cart(3,1:drawTimeIter:end);
text(xDraw,yDraw,zDraw,tLabelCell)
scatter3(xDraw,yDraw,zDraw,20,'b','*');
hold off;
%}

%% Test learning with basis functions
%{
example = 1;
[W,phi,q_reg,qd_reg,qdd_reg] = regress_basis_fncs(5,...
                                t_train{example},q_train{example},10);
figure('Name','Regressing on basis functs');
for i = 1:7
    subplot(7,3,3*(i-1)+1);
    plot(t_train{example},q_reg(:,i));
    hold on;
    plot(t_train{example},q_train{1}(:,i));
    subplot(7,3,3*(i-1)+2);
    plot(t_train{example},qd_reg(:,i));
    subplot(7,3,3*(i-1)+3);
    plot(t_train{example},qdd_reg(:,i));
end
% report maximum acc
loss_basis = norm(q_train{example} - q_reg,'fro')
max_acc_basis = max(max(abs(qdd_reg)))
%}

%% Train DMPs and evolve one

example = 15;
n_bf = 10;
tau = 1.0; % speed up/slow down movement
safe = false; %Jens' modification for safe-acc at the start of movement
cutoff = 10; % cutoff frequency
dmps = trainMultiDMPs(safe,cutoff,t_train{example},q_train{example},'d',n_bf);
dmps(1).can.tau = tau;
t_evolve = 1.0/tau;
N_evolve = t_evolve/dt;
q_dmp = zeros(7,N_evolve);
qd_dmp = zeros(7,N_evolve);
qdd_dmp = zeros(7,N_evolve);
ref = zeros(2*7,N_evolve);
for i = 1:7
    dmps(i).y0 = [q_train{example}(1,i);0;0];
    [x,y] = dmps(i).evolve(N_evolve);
    dmps(i).resetStates();
    q_dmp(i,:) = y(1,:);
    qd_dmp(i,:) = y(2,:);
    qdd_dmp(i,:) = y(3,:);
    ref(i,:) = y(1,:)';
    ref(i+7,:) = y(2,:)';
end
loss_dmp = norm(q_train{1} - q_dmp(:,1:end-1)','fro')
max_acc_dmp = max(max(abs(qdd_dmp)))
t_dmp = dt*(1:N_evolve);
figure('Name','DMP evolve');
for i = 1:7
    subplot(7,3,3*(i-1)+1);
    plot(t_dmp,q_dmp(i,:));
    %%{
    hold on;
    len = size(q_train{example},1);
    plot(t_dmp(1:len),q_train{example}(:,i));
    %legend('dmp','training');
    hold off;
    %}
    subplot(7,3,3*(i-1)+2);
    plot(t_dmp,qd_dmp(i,:));
    subplot(7,3,3*(i-1)+3);
    plot(t_dmp,qdd_dmp(i,:));
end

%% Visualize DMP in Cartesian space

wam = init_wam();
qd_plot = diff(q_train{example})/dt;
qd_plot = [qd_plot; qd_plot(end,:)];
ref_plot = [q_train{example}'; qd_plot'];
train_cart = wam.kinematics(ref_plot);

% draw Cartesian output
dmp_cart = wam.kinematics(ref);
figure('Name','DMP cartesian');
hold on;
grid on;
plot3(train_cart(1,:),train_cart(2,:),train_cart(3,:),'k-');
plot3(dmp_cart(1,:),dmp_cart(2,:),dmp_cart(3,:),'r--');
legend('training','dmp');
xlabel('x'); ylabel('y'); zlabel('z');
drawTimeIter = 40;
tLabel = t(1:drawTimeIter:N_evolve);
precision = 4;
tLabelCell = num2cell(tLabel,precision);
for i = 1:length(tLabelCell)
    tLabelCell{i} = num2str(tLabelCell{i});
end
% annotate some of the ball positions
xDraw = dmp_cart(1,1:drawTimeIter:end);
yDraw = dmp_cart(2,1:drawTimeIter:end);
zDraw = dmp_cart(3,1:drawTimeIter:end);
text(xDraw,yDraw,zDraw,tLabelCell)
scatter3(xDraw,yDraw,zDraw,20,'b','*');
hold off;

%% Compare safe dmps with unsafe dmps
%{
n_bf = 10;
tau = 1.0; % speed up/slow down movement
safe = true; %Jens' modification for safe-acc at the start of movement
dmps_safe = trainMultiDMPs(safe,t_train,q_train,'d',n_bf);
dmps_safe(1).can.tau = tau;
t_evolve = 1.0/tau;
N_evolve = t_evolve/dt;
q_dmp = zeros(7,N_evolve);
qd_dmp = zeros(7,N_evolve);
qdd_dmp = zeros(7,N_evolve);
ref_safe = zeros(2*7,N_evolve);
for i = 1:7
    dmps_safe(i).y0 = [q_train{1}(1,i);0;0];
    [x,y] = dmps_safe(i).evolve(N_evolve);
    dmps_safe(i).resetStates();
    q_dmp(i,:) = y(1,:);
    qd_dmp(i,:) = y(2,:);
    qdd_dmp(i,:) = y(3,:);
    ref_safe(i,:) = y(1,:)';
    ref_safe(i+7,:) = y(2,:)';
end
t_dmp = dt*(1:N_evolve);
figure('Name','DMP evolve');
for i = 1:7
    subplot(7,3,3*(i-1)+1);
    plot(t_dmp,q_dmp(i,:));
    subplot(7,3,3*(i-1)+2);
    plot(t_dmp,qd_dmp(i,:));
    subplot(7,3,3*(i-1)+3);
    plot(t_dmp,qdd_dmp(i,:));
end
% Visualize DMP in Cartesian space

% draw Cartesian output
dmp_cart_safe = wam.kinematics(ref_safe);
figure('Name','Save cart vs. unsafe cart');
hold on;
grid on;
plot3(train_cart(1,:),train_cart(2,:),train_cart(3,:),'r--');
plot3(dmp_cart(1,:),dmp_cart(2,:),dmp_cart(3,:),'b');
plot3(dmp_cart_safe(1,:),dmp_cart_safe(2,:),dmp_cart_safe(3,:),'k');
hold off;
%}

%% Save DMPs in JSON format
%%{
dmp_save = struct();
dmp_save.tau = 1.0;
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
fid = fopen(['../json/dmp_', int2str(example), '_', date, '.json'],'wt');
fwrite(fid, json_txt, 'char');
fclose(fid);
%}