% 复现网络重构  直接用差分法复现
% 结果: (1)ode45所使用的差分位置不确定->(2)可以用线性插值进行修复但插值的误差太大
% 希望: 自己写龙格库塔, 结合真实轨迹, 得到中心差分和导数, 在已知数值导数的位置上, 求解数值积分->main_9.m

clear
clc

% y0 = zeros(1,29); 
Y = ones(1,29); 

% (1)验证ode45可以固定步长来求解
step_length = 0.01; % 设置步长
[t,y] = ode45('fun_1_lorenz_solver',[0:step_length:2000],Y); % 求解 - 固定步长

% (2)计算x的数值导数, 基于差分法
[t,real_tr] = ode45('fun_0_lorenz',[0:0.00001:20],[10,8/3,28]);
dot_x = diff(real_tr(:,1))/0.0001;
save('mat_6_dot_x.mat', 'dot_x');

% (3)把差分结果传入ode45的求解函数中, 根据时间进行
step_length = 0.01; % 设置步长
[t,y] = ode45('fun_6_5_lorenz_solver',[0:step_length:0.1],Y); % 求解 - 固定步长

% 输出参数的结果
y(end,[25,26,27])

% 画出参数的运动
plot(y(:,25)) % a
hold on
plot(y(:,26)) % b
hold on 
plot(y(:,27)) % r
hold on 
legend('a','b','r') 