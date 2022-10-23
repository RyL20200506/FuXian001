% 尝试画正反lorenz的图像
close all
clear
clc


%% 还原失败
y0 = [1,1,1]; 
% 龙格库塔: 直接调用
length = 2;
[t,y] = fun_8_RungeKutta(@fun_0_lorenz,0,0.01,length,y0);
y = y';

% 画图
figure
plot(y(:,1),y(:,3));  % xz相图轨迹
grid on;


% 反向 龙格库塔 调用函数
% 最后一个点
the_point = y(end, :);
[t_n, y_n] = fun_8_RungeKutta(@fun_18_negative_lorenz, 0, 0.01, length, the_point);
y_n = y_n';

figure
back_n = length/0.01
plot(y_n(1:back_n ,1), y_n(1:back_n ,3));  % xz相图轨迹
grid on;

stop
%% 成功还原
close all
clear
clc

y0 = [1,1,1]; 
length = 200;

% 龙格库塔: 函数展开
derivatives = [];  % 用于记录导数
ufunc = @fun_0_lorenz;
a = 0;
h = 0.01;
b = length;
n=floor((b-a)/h);       %步数
x(1)=a;                 %时间起点
y(:,1)=y0;              %赋初值，可以是向量，但是要注意维数
for i=1:n               %龙格库塔方法进行数值求解
    x(i+1)=x(i)+h;
    k1=ufunc(x(i),y(:,i));
    k2=ufunc(x(i)+h/2,y(:,i)+h*k1/2);
    k3=ufunc(x(i)+h/2,y(:,i)+h*k2/2);
    k4=ufunc(x(i)+h,y(:,i)+h*k3);
    y(:,i+1)=y(:,i)+h*(k1+2*k2+2*k3+k4)/6;
    derivatives = [derivatives, (k1+2*k2+2*k3+k4)/6];  % 将导数的结果进行记录.
end
y = y';


% 画图
figure
plot(y(:,1),y(:,3));  % xz相图轨迹
grid on;

% -----------------------------

% 最后一个点
the_point = y(end, :);

% 反向 龙格库塔 调用函数
[t_n, y_n] = fun_8_RungeKutta(@fun_18_negative_lorenz, 0, 0.01, length, the_point);
y_n = y_n';

% 反向 基于导数反推
y_n = [the_point];
y_n = y_n';  % 转为时间在列
for i=1:size(derivatives, 2)
    y_n(:, i+1) = y_n(:,i)-h*derivatives(:,end-(i-1));  % (1)要将导数方向修改 (2)从最后一个的导数开始执行
end
y_n = y_n';  % 转化为时间在行

figure
back_n = size(derivatives, 2);
plot(y_n(1:back_n ,1), y_n(1:back_n ,3));  % xz相图轨迹
grid on;

% start point
start_point = y_n(end, :);