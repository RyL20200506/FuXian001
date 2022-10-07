% 复现网络重构
clear
clc

% % 获取真实轨迹 - 法1: 基于ode45
% Y = [-8, 7, 27, ones(1,26)]; 
% [t,y] = ode45('fun_1_lorenz_solver',[0,2000],Y); % 求解

% 获取真实轨迹 - 法2: 自定义龙格库塔
Y = [-8, 7, 27, ones(1,26)]; 
end_point = 5000;
[t,y] = fun_8_RungeKutta(@fun_1_lorenz_solver,0,0.005,end_point,Y);
X_n = y';

% 画出参数的运动
figure
plot(X_n(:,25)) % a
hold on 
plot(X_n(:,26)) % b
hold on 
plot(X_n(:,27)) % r
hold on 
legend('a','b','r') 
% save('mat_1_y', 'y')
para = X_n(end, [25,26,27])



%% 
Y = [-8, 7, 27, ones(1,28)]; 
end_point = 5000;

% 法2
[t,X_n] = ode45('fun_1_lorenz_solver1',[0,5000], Y); % 求解
para = X_n(end, [27,28,29])

% 画出参数的运动
figure
plot(X_n(:,27)) % a
hold on 
plot(X_n(:,28)) % b
hold on 
plot(X_n(:,29)) % r
hold on 
legend('a','b','r') 
% save('mat_1_y', 'y')

% 预测
para = X_n(end, [27,28,29])
start_to_pred_point = X_n(end,1:3)  % 结果: 大概能预测100个点
% start_to_pred_point = X_n(1,1:3)  % 结论: 大概能预测500个点, 能够预测的时间长度, 与从什么时候开始预测也有关
[t,y] = fun_8_RungeKutta(@fun_0_lorenz,0,0.01,end_point, start_to_pred_point);
X_n_1 = y';
[t,y] = fun_8_RungeKutta_para(@fun_0_lorenz_para, 0, 0.01, end_point, start_to_pred_point, para);
X_pred = y';

% 画出相图轨迹
figure
plot(X_n_1(:,2)', 'linewidth', 1.5);
hold on
plot(X_pred(:,2),'--', 'linewidth',1.5);
xlim([0 3000]);

h1 = legend('real','predict') ;
xlabel('\it Step \rm', 'fontsize',14);
ylabel('\it Value \rm', 'fontsize',14);
% 自定义刻度
xtickformat('%.1f');
ytickformat('%.1f');
ax = gca;
% ax.XAxis.Exponent = 5;
% ylim([0 30])
% 调整字体大小
set(gca,'FontSize',14)  %是设置刻度字体大小
set (h1,'box','off')

