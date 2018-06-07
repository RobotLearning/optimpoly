%% Train multi-dmps
% pattern determines rhythmic or discrete pattern generation
function dmp = trainMultiDMPs(t,q,pat,n_bf)

f = 500; %200 Hz recording

% number of demonstrations
D = length(q);
dof = size(q{1},2);
tf = t{end};

dt = tf(2) - tf(1);
t_common = linspace(dt,tf(end)-tf(1),length(tf)); 
goals = zeros(D,dof);
inits = zeros(D,dof);

for i = 1:D
    qd{i} = diff(q{i})/dt;
    t{i} = t{i} - t{i}(1);
    td = t{i}(1:end-1);
    qdd{i} = diff(qd{i})/dt;
    tdd = td(1:end-1);
    q{i} = interp1(t{i},q{i},t_common,'linear','extrap');
    qd{i} = interp1(td,qd{i},t_common,'linear','extrap');
    qdd{i} = interp1(tdd,qdd{i},t_common,'linear','extrap');
end
    

for i = 1:D
    
    if strcmp(pat,'d')
        goals(i,:) = q{i}(end,:);
        vels(i,:) = qd{i}(end,:);
        accs(i,:) = qdd{i}(end,:);
    else
        goals(i,:) = (min(q{i}) + max(q{i})) / 2;
    end
    inits(i,:) = q{i}(1,:);
end

% canonical system
h = 1/f; % 200 Hz recordings
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
    q_reg = [q_reg; q{i}];
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
        dmp(i) = DDMP(can,alpha,beta,g(i),yin(:,i));
    else
        amp(i) = 1;
        dmp(i) = RDMP(can,alpha,beta,g(i),amp(i),yin(:,i));
    end
    dmp(i).regressLive(q_reg(:,i),qd_reg(:,i),qdd_reg(:,i),goals(:,i));
    
end

disp('Initializing DMP to yin = ');
disp(yin(1,:)');
disp('Goal state for DMP, goal = ');
disp(g(:));