% 复现: Wang - 重构动力学
% Wang 误差函数

close all;clear;clc;

% 构造训练集
y0 = [1,1,1];  % 初始值
[t, real_trajectory] = ode45('fun_0_lorenz',[0,100],y0); % 带入之前的函数求解真实轨迹
X = real_trajectory(1:end-1, :);  % 分割出特征
Y = real_trajectory(2:end, :);  % 再分割出待预测标签

% 构造需要最小化的函数, x是希望拟合的参数, 输入tdata进去, 希望得到ydata
fun = @(PARA)fun_3_costfun(PARA,X,Y);

% 搜索最优参数
PARA_0 = rand(48,1);  % 初始参数
options = optimset('MaxFunEvals',10000000);
best_PARA = fminsearch(fun,PARA_0, options);
% 更新迭代了几次
best_PARA = fminsearch(fun,best_PARA, options);  
best_PARA = fminsearch(fun,best_PARA, options); 
best_PARA = fminsearch(fun,best_PARA, options); 
best_PARA = fminsearch(fun,best_PARA, options); 

% --------------------------------------------------
% 更换训练集
% 构造训练集
y0 = [1,1,1];  % 初始值
[t, real_trajectory] = ode45('fun_0_lorenz',[20,50],y0); % 带入之前的函数求解真实轨迹
X = real_trajectory(1:end-1, :);  % 分割出特征
Y = real_trajectory(2:end, :);  % 再分割出待预测标签

% 构造需要最小化的函数, x是希望拟合的参数, 输入tdata进去, 希望得到ydata
fun = @(PARA)fun_3_costfun(PARA,X,Y);

% 搜索最优参数
PARA_0 = rand(48,1);  % 初始参数
options = optimset('MaxFunEvals',10000000);
% 更新迭代几次
best_PARA = fminsearch(fun,best_PARA, options);  
best_PARA = fminsearch(fun,best_PARA, options); 


% --------------------------------------------------
% --------------------------------------------------
% 检查拟合质量: 就看看模型的精度 
X0 = X(1,:);
predict_result = fun_3_4_predict(best_PARA, X0, 1000); 

% 
plot(predict_result(:, 1), '*');
hold on
plot(X(1:1000, 1),'r');

% 看看best para的拟合效果
A = bestx(1);
lambda = bestx(2);
yfit = A*exp(-lambda*tdata);
plot(tdata,ydata,'*');
hold on
plot(tdata,yfit,'r');
xlabel('tdata')
ylabel('Response Data and Curve')
title('Data and Best Fitting Exponential Curve')
legend('Data','Fitted Curve')
hold off
