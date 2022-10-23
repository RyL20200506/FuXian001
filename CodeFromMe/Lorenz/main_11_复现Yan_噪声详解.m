% 复现Yan. 在main_10基础上将噪声修改为, 每一步的噪声都会影响到下一时刻的结果值.
clear
clc

step = 0.01;  % 拟合步长为step, 真实轨迹的步长为step/2  (因为:RungeKutta会用到1/2的步长去计算)
real_b = 20;  % 真实轨迹的结束时间 
real_a = 0;  % 真实轨迹的起始时间 
y0 = [-8,7,27];  % 真实轨迹的初始点
interp_step = step  % step 或 step/2: #注: 不可以用step/2, 会有很严重的问题
interp_a = real_a + interp_step;  % 插值F的开始时间  取step/2的中心差分后, 前面会少1个时刻的数据;
interp_b = real_b - interp_step;  % 插值F的结束时间
D = 0.2;

% 数据准备: 获取已知可观察轨迹, 但是是有噪声版本
a=real_a;
b=real_b;
h=interp_step;
n=floor((b-a)/h);       %步数
x(1)=a;                 %时间起点
y(:,1)=y0;              %赋初值，可以是向量，但是要注意维数
rng(1);                 % 设置随机数的种子
for i=1:n               %龙格库塔方法进行数值求解
    x(i+1)=x(i)+h;
    k1=fun_0_lorenz(x(i),y(:,i));
    k2=fun_0_lorenz(x(i)+h/2,y(:,i)+h*k1/2);
    k3=fun_0_lorenz(x(i)+h/2,y(:,i)+h*k2/2);
    k4=fun_0_lorenz(x(i)+h,y(:,i)+h*k3);
    y(:,i+1)=y(:,i)+h*(k1+2*k2+2*k3+k4)/6 + D*randn(1);
end
real_trajectory = y';

% 画出原始的真实轨迹图
figure
plot3(real_trajectory(:,1),real_trajectory(:,2),real_trajectory(:,3));  % 画图
grid on;  % 画格子

% 只画x-z轴的图像
plot(real_trajectory(:,1),real_trajectory(:,3));  % xz相图轨迹

% 注: 苛刻的说, 上面画的图象, 并不是论文中希望画出的图像; 主要原因在于, 原文是在dot{x}上加噪声, 而我们是在结果上加噪声
% 这要求我们最好重新写一个带噪声的Lorenz函数

% 用公式(1)去画噪声.
[t,y] = fun_8_RungeKutta(@fun_11_lorenz_noise,1,interp_step,100,y0); 
real_trajectory = y';
figure
plot(real_trajectory(:,1),real_trajectory(:,3));  % xz相图轨迹
