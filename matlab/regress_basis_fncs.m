%% Regress movement using basis functions for each dof
function [W,phi,Q,Qd,Qdd] = regress_basis_fncs(filt_freq,t,q,n_bf)

fc = filt_freq;
fs = 500; % Hz recording

[b,a] = myButter2ndOrder(fc/(fs/2));
q = filtfilt(b,a,q);
t = t - t(1);
phi = form_phi(t,n_bf); % add intercept and linear trend
W = zeros(size(phi,2),7);
lambda = 1;
for i = 1:7
    %W(:,i) = (phi'*phi + lambda*eye(n_bf+2))\(phi'*q{1}(:,i));
    W(:,i) = phi \ q(:,i);
end
dt = t(2)-t(1); %assumed uniform t
Q = phi*W;
Qd = diff(Q)/dt;
Qdd = diff(Qd)/dt;
Qd = interp1(t(1:end-1),Qd,t);
Qdd = interp1(t(1:end-2),Qdd,t);

end

function phi = form_phi(t,nbf)

N = length(t);
%phi = zeros(N,nbf);
phi = zeros(N,nbf+2);
c = linspace(t(end)/nbf,t(end),nbf);
h = 5*ones(1,nbf);
for i = 1:N
    %phi(i,:) = basis_fncs(t(i),c,h);
    phi(i,:) = [1,t(i),basis_fncs(t(i),c,h)];
end

end

function y = basis_fncs(x,c,h)

    y = exp(-h.*((x-c).^2));
    
end