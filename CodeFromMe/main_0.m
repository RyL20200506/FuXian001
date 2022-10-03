% 用于复现网络重构
clear
clc

y0 = [1,1,1]; 
[t,y] = ode45('fun_0_lorenz',[0,200],y0); % 求解
% 画图
plot(y(:,1),y(:,3));  % xz相图轨迹
grid on;
