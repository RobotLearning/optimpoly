%% Train multi-dmps
% pattern determines rhythmic or discrete pattern generation
function dmp = trainMultiDMPs(safe,filt_freq,t,q,pat,n_bf)

fc = filt_freq;
fs = 500; % Hz recording

% number of demonstrations
if iscell(q)
    D = length(qs);
    dof = size(qs{1},2);
    ts = t;
    qs = q;
else
    D = 1;
    dof = size(q,2);
    ts{1} = t;
    qs{1} = q;
end
    
tf = ts{end};

dt = tf(2) - tf(1);
t_common = linspace(dt,tf(end)-tf(1),length(tf)); 
goals = zeros(D,dof);
inits = zeros(D,dof);
[b,a] = myButter2ndOrder(fc/(fs/2));
%filt_order = 2;
%[b,a] = butter(filt_order,fc/(fs/2));

for i = 1:D
    qs{i} = filtfilt(b,a,qs{i});
    qd{i} = diff(qs{i})/dt;
    ts{i} = ts{i} - ts{i}(1);
    td = ts{i}(1:end-1);
    qdd{i} = diff(qd{i})/dt;
    tdd = td(1:end-1);
    qs{i} = interp1(ts{i},qs{i},t_common,'linear','extrap');
    qd{i} = interp1(td,qd{i},t_common,'linear','extrap');
    qdd{i} = interp1(tdd,qdd{i},t_common,'linear','extrap');
end
    

for i = 1:D
    
    if strcmp(pat,'d')
        goals(i,:) = qs{i}(end,:);
        vels(i,:) = qd{i}(end,:);
        accs(i,:) = qdd{i}(end,:);
    else
        goals(i,:) = (min(qs{i}) + max(qs{i})) / 2;
    end
    inits(i,:) = qs{i}(1,:);
end

% canonical system
h = 1/fs; 
tau = 1/t_common(end);
alpha = 25;
beta = alpha/4;
ax = 1;
% number of basis functions
numbf = n_bf;
can = CAN(h,ax,tau,numbf,t_common(end),pat);

% stack cells 
q_reg = [];
qd_reg = [];
qdd_reg = [];
for i = 1:D
    q_reg = [q_reg; qs{i}];
    qd_reg = [qd_reg; qd{i}];
    qdd_reg = [qdd_reg; qdd{i}];
end
g = zeros(1,dof);
yin = zeros(3,dof);

for i = 1:dof
    
    % tau is fixed to be one
    % goal and amplitude are initialized here
    % these are not important as dmps can be extended to any yin and goal
    yin(:,i) = [inits(1,i);0;0];
    g(i) = goals(1,i);
    % initial states of DMPs
    if strcmp(pat,'d')
        dmp(i) = DDMP(safe,can,alpha,beta,g(i),yin(:,i));
    else
        amp(i) = 1;
        dmp(i) = RDMP(safe,can,alpha,beta,g(i),amp(i),yin(:,i));
    end
    dmp(i).regressLive(q_reg(:,i),qd_reg(:,i),qdd_reg(:,i),goals(:,i));
    
end

disp('Initializing DMP to yin = ');
disp(yin(1,:)');
disp('Goal state for DMP, goal = ');
disp(g(:));