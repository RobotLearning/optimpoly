% Plot joints act,des data saved in SL

clc; clear; close all;
M = dlmread('../joints.txt');
if (size(M,2) == 14)
    % in this case qdes was also saved
    plot_qdes = true;
end

t1 = 0.002:0.002:1.0;
figure;
for i = 1:7
    subplot(7,1,i); 
    plot(t1,M(1:500,i));
    hold on;
    if plot_qdes
        plot(t1,M(1:500,i+7));
    end
end
N2 = size(M,1) - 500;
t2 = 0.002 * (1:N2);
figure;
for i = 1:7
    subplot(7,1,i); 
    plot(t2,M(501:end,i)); 
    hold on;
    if plot_qdes
        plot(t2,M(501:end,i+7));
    end
end