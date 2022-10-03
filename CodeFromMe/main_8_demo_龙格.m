% ͨ����������ķ���, ����ode45
% �ο�: https://www.cnblogs.com/ybqjymy/p/13645478.html
clear all;
close all;
clc;

%ϵͳ���������
[t,h] = ode45('fun_0_lorenz',[0:0.01: 40],[12 4 0]);  % ���뺯��, ��Χ, ��ֵ
plot3(h(:,1),h(:,2),h(:,3));
grid on;

%�Զ������������
[t1,h1]=fun_8_RungeKutta(@fun_0_lorenz,0, 0.01,40, [12 4 0]);
figure;
plot3(h1(1,:),h1(2,:),h1(3,:),'r')
grid on;


