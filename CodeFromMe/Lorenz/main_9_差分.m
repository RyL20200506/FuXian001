% 复现Yan. 使用自己实现的龙格库塔, 去结合差分法, 进行求解
clear
clc

step = 0.01;  % 拟合步长为step, 真实轨迹的步长为step/2  (因为:RungeKutta会用到1/2的步长去计算)
real_b = 2000;  % 真实轨迹的结束时间 
real_a = 0;  % 真实轨迹的起始时间 
y0 = [-8,7,27];  % 真实轨迹的初始点 
interp_step = step  % step 或 step/2: #注: 不可以用step/2, 会有很严重的问题
interp_a = real_a + interp_step;  % 插值F的开始时间  取step/2的中心差分后, 前面会少1个时刻的数据;
interp_b = real_b - interp_step;  % 插值F的结束时间

% 数据准备: 获取已知可观察轨迹
[t,real_trajectory] = fun_8_RungeKutta(@fun_0_lorenz,real_a,interp_step,real_b,y0);  % 求解轨迹: 固定步长: 步长是拟合所用步长的一半
real_trajectory = real_trajectory';
known_X = real_trajectory(:, 1);  % 真实且已知的轨迹
% 插值: 直接插X的值
FX = griddedInterpolant(t, known_X);

% 数据准备: 获取已知可观察轨迹的中心差分, 从而
diff_before_X = known_X(2:end-1) - known_X(1:end-2);  % 前向差分: #验证无误
diff_after_X = (known_X(3:end) - known_X(2:end-1));  % 后向差分  #验证无误
diff_center_X = (diff_before_X+diff_after_X)/2;  % 中心差分
dot_X = diff_center_X/(interp_step);  % 数值导数
save('mat_9_dot_X', 'dot_X');  % 保存中间变量

% 定义插值函数 fun_9_lorenz_solver: 使其能接受外部输入dot_X数据, 并能根据step选取特定的dot_X的数据
time_range = [interp_a: interp_step: interp_b];
F = griddedInterpolant(time_range, dot_X);

% 开始求解
n=floor((real_b-real_a)/step);       % 步数: 只要步数能够被dot_X覆盖, 就不会有误差
time(1)=0;             % 时间列表
y=[-8; 7; 27; ones(26,1)];              % 赋初值，可以是向量，但是要注意维数

% n=1000
for i=1:n          % 龙格库塔方法进行数值求解
    time(i+1)=time(i)+step;
    k1=fun_9_lorenz_solver(time(i),y(:,i),F, FX);
    k2=fun_9_lorenz_solver(time(i)+step/2,y(:,i)+step*k1/2,F, FX);
    k3=fun_9_lorenz_solver(time(i)+step/2,y(:,i)+step*k2/2,F, FX);
    k4=fun_9_lorenz_solver(time(i)+step,y(:,i)+step*k3,F, FX);
    y(:,i+1)=y(:,i)+step*(k1+2*k2+2*k3+k4)/6;
end
result = y';
result(end,[25,26,27])  % 这是结果


% 画出参数的运动图
figure
plot(result(:,25)) % a
hold on 
plot(result(:,26)) % b
hold on 
plot(result(:,27)) % r
hold on 
legend('a','b','r') 

stop  % 下面都是调试

% 画出拟合步骤的轨迹图
figure;
plot3(result(:,1),result(:,2),result(:,3),'r')
grid on;

% 画出原始的真实轨迹图
figure
plot3(real_trajectory(:,1),real_trajectory(:,2),real_trajectory(:,3));  % 画图
grid on;  % 画格子

% 保存轨迹结果
save('mat_9_result','result')

% 目前的困难: F的值(3)与微分方程的真实导数不一样, 且F的值与dot_X是对应的.
% (1)调试: 对比两种不同真实轨迹, 3个变量, 和29个变量的轨迹导数是否一致-> 下面代码表明它们的结果完全一致
[t,real_3] = fun_8_RungeKutta(@fun_9_lorenz,real_a,step/2,real_b,[-8,7,27]);  % 求解轨迹: 3个元素的轨迹
real_3 = real_3';
known_X_3 = real_3(:, 1);  % 真实且已知的轨迹
[t,real_29] = fun_8_RungeKutta(@fun_9_lorenz_solver,real_a,step/2,real_b,[-8; 7; 27; ones(26,1)]);  % 求解轨迹: 29个元素的轨迹
real_29 = real_29';
known_X_29 = real_29(:, 1);  % 与上面完全一致.

% 调试: 对比差分法得到的导数, 与真实轨迹的导数是否一致: 感觉可能不一致
fun_8_RungeKutta(@fun_9_lorenz,real_a,step/2,real_b,y0);

% (2)调试: 在步长不相同的情况下, 对dotX结果的影响; 因为真实数据的步长是step/2, 虚拟数据的步长是step 
% -> 确实是这个问题导致的! 步长会对dotX造成很大的影响, 它们的轨迹应该都有很大区别!
% -> Q:那么, dot_X应该如何得到呢?
[t,real_3] = fun_8_RungeKutta(@fun_9_lorenz,real_a, step/2, real_b,[-8,7,27]);  % 求解轨迹: 3个元素的轨迹
real_3 = real_3';
known_X_3 = real_3(:, 1);  % 真实且已知的轨迹
diff_center_X_3 = fun_9_diff_center(known_X_3);
dot_X_3 = diff_center_X_3/(step/2);
time_range = [real_a+step/2: step/2: real_b-step/2];
F_3 = griddedInterpolant(time_range, dot_X_3);

[t,real_29] = fun_8_RungeKutta(@fun_9_lorenz_solver,real_a, step, real_b,[-8; 7; 27; ones(26,1)]);  % 求解轨迹: 29个元素的轨迹
real_29 = real_29';
known_X_29 = real_29(:, 1);  % 真实且已知的轨迹
diff_center_X_29 = fun_9_diff_center(known_X_29);
dot_X_29 = diff_center_X_29/(step);
time_range = [real_a+step: step: real_b-step];
F_29 = griddedInterpolant(time_range, dot_X_29);

% (3)调试: 由此可知, 我们将用相同步长, 再用线性插值的方法去完成
% -> 