% ���ڸ��������ع�
clear
clc

y0 = [1,1,1]; 
[t,y] = ode45('fun_0_lorenz',[0,200],y0); % ���
% ��ͼ
plot(y(:,1),y(:,3));  % xz��ͼ�켣
grid on;
