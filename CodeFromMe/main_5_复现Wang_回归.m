% 复现: Wang - 重构动力学: 用线性回归来求解

close all;clear;clc;

% 构造已知集合
y0 = [-8,7,27];  % 初始值
[t, real_trajectory] = ode45('fun_0_lorenz', [0,100], y0);  % 带入之前的函数求解真实轨迹

% 一些参数
t_unit = 0.01; 

% 准备数据集
% (1)取出变量
X = real_trajectory(:,1);  % 取出列向量, 时间随行数增加
Y = real_trajectory(:,2);
Z = real_trajectory(:,3);

% (2)准备: 可能需要的基向量
XY = X.*Y; 
XZ = X.*Z; 
YZ = Y.*Z; 
XYZ = XY.*Z;

% (3)合并
about_all_x = [ones(size(X)), X, X*t_unit, Y*t_unit, Z*t_unit, XY*t_unit, XZ*t_unit, YZ*t_unit, XYZ*t_unit, X*(t_unit^2), Y*(t_unit^2), Z*(t_unit^2), XY*(t_unit^2), XZ*(t_unit^2), YZ*(t_unit^2), XYZ*(t_unit^2)];
about_all_y = [ones(size(Y)), Y, X*t_unit, Y*t_unit, Z*t_unit, XY*t_unit, XZ*t_unit, YZ*t_unit, XYZ*t_unit, X*(t_unit^2), Y*(t_unit^2), Z*(t_unit^2), XY*(t_unit^2), XZ*(t_unit^2), YZ*(t_unit^2), XYZ*(t_unit^2)];
about_all_z = [ones(size(Z)), Z, X*t_unit, Y*t_unit, Z*t_unit, XY*t_unit, XZ*t_unit, YZ*t_unit, XYZ*t_unit, X*(t_unit^2), Y*(t_unit^2), Z*(t_unit^2), XY*(t_unit^2), XZ*(t_unit^2), YZ*(t_unit^2), XYZ*(t_unit^2)];

% 训练
known_all_x = about_all_x(1:end-1,:);  % 因为最后一行只充当label
need_all_x = X(2:end);  % 因为第一行只充当feature
para_x = regress(need_all_x, known_all_x)  % 线性回归

known_all_y = about_all_y(1:end-1,:);  % 因为最后一行只充当label
need_all_y = Y(2:end);  % 因为第一行只充当feature
para_y = regress(need_all_y, known_all_y)  % 线性回归

known_all_z = about_all_z(1:end-1,:);  % 因为最后一行只充当label
need_all_z = Z(2:end);  % 因为第一行只充当feature
para_z = regress(need_all_z, known_all_z)  % 线性回归

% -----------------------------------------

% 查看训练效果
% 迭代生成X的结果
the_known_x = known_all_x(1,:);
the_knwon_y = known_all_y(1,:);
the_known_z = known_all_z(1,:);

% 
X_predict_history = [];  % 容器 记录结果
Y_predict_history = []; 
Z_predict_history = []; 

for i=1:size(known_all_x)
    x = the_known_x*para_x;  % 
    y = the_knwon_y*para_y;
    z = the_known_z*para_z;
    
    X_predict_history = [X_predict_history; x];  % 列合并 
    Y_predict_history = [Y_predict_history; y]; 
    Z_predict_history = [Z_predict_history; z]; 
    
    % 更新基向量
    the_known_x = [1, x, x*t_unit, y*t_unit, z*t_unit, (x*y)*t_unit, (x*z)*t_unit, (y*z)*t_unit, (x*y*z)*t_unit, x*(t_unit^2), y*(t_unit^2), z*(t_unit^2), x*y*(t_unit^2), x*z*(t_unit^2), y*z*(t_unit^2), x*y*z*(t_unit^2)];
    the_known_y = [1, y, x*t_unit, y*t_unit, z*t_unit, (x*y)*t_unit, (x*z)*t_unit, (y*z)*t_unit, (x*y*z)*t_unit, x*(t_unit^2), y*(t_unit^2), z*(t_unit^2), x*y*(t_unit^2), x*z*(t_unit^2), y*z*(t_unit^2), x*y*z*(t_unit^2)];
    the_known_z = [1, z, x*t_unit, y*t_unit, z*t_unit, (x*y)*t_unit, (x*z)*t_unit, (y*z)*t_unit, (x*y*z)*t_unit, x*(t_unit^2), y*(t_unit^2), z*(t_unit^2), x*y*(t_unit^2), x*z*(t_unit^2), y*z*(t_unit^2), x*y*z*(t_unit^2)];
end

% 画图
plot(X_predict_history, 'r');
hold on

% y0 = [1,1,1];  % 初始值
% [t, real_trajectory] = ode45('fun_0_lorenz', [0,200], y0);  % 带入之前的函数求解真实轨迹
% plot(real_trajectory(:,1), '*');
% hold on
% 
% % -----------------------------------------------
% % 貌似很好, 其实是错误的结果
% the_known_x = known_all_x(1,:);
% the_knwon_y = known_all_y(1,:);
% the_known_z = known_all_z(1,:);
% 
% % 
% X_predict_history = [];  % 容器 记录结果
% Y_predict_history = []; 
% Z_predict_history = []; 
% 
% for i=1:size(known_all_x)
%     x = the_known_x*para_x;  % 
%     y = the_knwon_y*para_y;  % 
%     z = the_known_z*para_z;  % 
%     
%     X_predict_history = [X_predict_history; x];  % 列合并 
%     Y_predict_history = [Y_predict_history; y];  % 
%     Z_predict_history = [Z_predict_history; z];  
%     
%     % 更新基向量
%     the_known_x = [1, x, x*t_unit, y*t_unit, z*t_unit, (x*y)*t_unit, (x*z)*t_unit, (y*z)*t_unit, (x*y*z)*t_unit, x*(t_unit^2), y*(t_unit^2), z*(t_unit^2), x*y*(t_unit^2), x*z*(t_unit^2), y*z*(t_unit^2), x*y*z*(t_unit^2)];
%     the_known_y = [1, y, x*t_unit, y*t_unit, z*t_unit, (x*y)*t_unit, (x*z)*t_unit, (y*z)*t_unit, (x*y*z)*t_unit, x*(t_unit^2), y*(t_unit^2), z*(t_unit^2), x*y*(t_unit^2), x*z*(t_unit^2), y*z*(t_unit^2), x*y*z*(t_unit^2)];
%     the_known_z = [1, z, x*t_unit, y*t_unit, z*t_unit, (x*y)*t_unit, (x*z)*t_unit, (y*z)*t_unit, (x*y*z)*t_unit, x*(t_unit^2), y*(t_unit^2), z*(t_unit^2), x*y*(t_unit^2), x*z*(t_unit^2), y*z*(t_unit^2), x*y*z*(t_unit^2)];
% end
% 
% % 画图
% plot(X_predict_history, 'r');
% hold on
% 
% 
