% 通过龙格库塔的方法, 代替ode45
% 参考: https://www.cnblogs.com/ybqjymy/p/13645478.html
clear all;
close all;
clc;

%系统龙格库塔法
[t,h] = ode45('fun_0_lorenz',[0:0.01: 40],[12 4 0]);  % 输入函数, 范围, 初值
plot3(h(:,1),h(:,2),h(:,3));
grid on;

%自定义龙格库塔法
[t1,h1]=fun_8_RungeKutta(@fun_0_lorenz,0, 0.01,40, [12 4 0]);
figure;
plot3(h1(1,:),h1(2,:),h1(3,:),'r')
grid on;


