%% Load data saved from python

clc; clear; close all;
save_tikz = false;
example = 0;
robot_file = ['python/move_robot_', num2str(example), '.txt'];
R = dlmread(robot_file);
N = size(R,1);
t_robot = R(:,1);
t_robot = t_robot - t_robot(1);
J = R(:,2:8);
X = R(:,9:11);

ball_file = 'python/move_ball_1.txt';
B = dlmread(ball_file);
t_ball = B(:,1);
t_ball = t_ball - t_ball(1);
B = B(:,2:end);

%% joint plot
figure;
for i = 1:7
     subplot(7,1,i);
     plot(t_robot,J(:,i),'LineWidth',1.0);
     ylabel(['q_', num2str(i)]);
end
xlabel('t');
%cleanfigure;
if save_tikz
    save_file = ['../learning-to-serve/Pictures/demo_kin_', num2str(example), '.tikz'];
    matlab2tikz(save_file,'height', '\fheight', 'width', '\fwidth');
    close all;
end
%% cartesian output

figure('Name','Training cartesian data');
hold on;
grid on;
scatter3(X(:,1),X(:,2),X(:,3),'r');
xlabel('x'); 
ylabel('y'); 
zlabel('z');
% add ball
draw_ball = 1:100;
scatter3(B(draw_ball,1),B(draw_ball,2),B(draw_ball,3),'b');
legend('robot','ball');
drawTimeIter = 40;
tLabel = t_robot(1:drawTimeIter:end);
precision = 4;
tLabelCell = num2cell(tLabel,precision);
for i = 1:length(tLabelCell)
    tLabelCell{i} = num2str(tLabelCell{i});
end
% annotate some of the ball positions
xDraw = X(1:drawTimeIter:end,1);
yDraw = X(1:drawTimeIter:end,2);
zDraw = X(1:drawTimeIter:end,3);
text(xDraw,yDraw,zDraw,tLabelCell)
scatter3(xDraw,yDraw,zDraw,4,'b','*');

tLabel = t_ball(1:drawTimeIter:draw_ball(end));
precision = 4;
tLabelCell = num2cell(tLabel,precision);
for i = 1:length(tLabelCell)
    tLabelCell{i} = num2str(tLabelCell{i});
end
xDraw = B(1:drawTimeIter:draw_ball(end),1);
yDraw = B(1:drawTimeIter:draw_ball(end),2);
zDraw = B(1:drawTimeIter:draw_ball(end),3);
text(xDraw,yDraw,zDraw,tLabelCell)
scatter3(xDraw,yDraw,zDraw,4,'r','*');

hold off;
%cleanfigure;
if save_tikz
    matlab2tikz('../learning-to-serve/Pictures/demo_kin2.tikz',...
        'height', '\fheight', 'width', '\fwidth');
    close all;
end