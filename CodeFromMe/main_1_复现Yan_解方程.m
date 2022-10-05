% ���������ع�
clear
clc

% % ��ȡ��ʵ�켣 - ��1: ����ode45
% Y = [-8, 7, 27, ones(1,26)]; 
% [t,y] = ode45('fun_1_lorenz_solver',[0,2000],Y); % ���

% ��ȡ��ʵ�켣 - ��2: �Զ����������
% Y = [-8, 7, 27, ones(1,26)]; 
% end_point = 5000;
% [t,y] = fun_8_RungeKutta(@fun_1_lorenz_solver1,0,0.005,end_point,Y);
% X_n = y';

Y = [-8, 7, 27, ones(1,28)]; 
end_point = 5000;

% ��1
% [t,y] = fun_8_RungeKutta(@fun_1_lorenz_solver1,0,0.01,end_point,Y);
% X_n = y';

% ��2
[t,X_n] = ode45('fun_1_lorenz_solver1',[0,5000], Y); % ���

% ������
X_n(end,27)
X_n(end,28)
X_n(end,29)

stop

% �����������˶�
figure
plot(X_n(:,25)) % a
hold on 
plot(X_n(:,26)) % b
hold on 
plot(X_n(:,27)) % r
hold on 
legend('a','b','r') 
save('mat_1_y', 'y')

% ��������Ľ��
para = X_n(end, [25,26,27])
% para = [10, 8/3, 28];
% para = [9.9979, 2.7158, 27.8544];

% �����������������ֵ���
start_to_pred_point = X_n(end,1:3)  % ���: �����Ԥ��100����
% start_to_pred_point = X_n(1,1:3)  % ����: �����Ԥ��500����, �ܹ�Ԥ���ʱ�䳤��, ���ʲôʱ��ʼԤ��Ҳ�й�
[t,y] = fun_8_RungeKutta(@fun_0_lorenz,0,0.01,end_point, start_to_pred_point);
X_n_1 = y';
[t,y] = fun_8_RungeKutta_para(@fun_0_lorenz_para, 0, 0.01, end_point, start_to_pred_point, para);
X_pred = y';

% ������ͼ�켣
figure
plot(X_n_1(:,2)', 'linewidth', 1);
hold on
plot(X_pred(:,2),'--', 'linewidth',1);
xlim([0 1000]);


