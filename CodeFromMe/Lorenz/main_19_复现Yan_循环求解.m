clear
clc
rng(1);
step = 0.01;  % 拟合步长为step, 真实轨迹的步长为step/2  (因为:RungeKutta会用到1/2的步长去计算)
real_b = 30;  % 真实轨迹的结束时间 
real_a = 0;  % 真实轨迹的起始时间 
y0 = [-8,7,27];  % 真实轨迹的初始点 
interp_step = step  % step 或 step/2: #注: 不可以用step/2, 会有很严重的问题
interp_a = real_a + interp_step;  % 插值F的开始时间  取step/2的中心差分后, 前面会少1个时刻的数据;
interp_b = real_b - interp_step;  % 插值F的结束时间

% 数据准备: 获取已知可观察轨迹
[t,y] = fun_8_RungeKutta(@fun_11_lorenz_noise, real_a, interp_step, real_b, y0); 
real_trajectory = y';
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
y = [-8; 7; 27; ones(28,1)];              % 赋初值，可以是向量，但是要注意维数
n=floor((real_b-real_a)/step);       % 步数: 只要步数能够被dot_X覆盖, 就不会有误差
time(1)=0;             % 时间列表
% n=1000
for i=1:n          % 龙格库塔方法进行数值求解
    time(i+1)=time(i)+step;
    k1=fun_17_lorenz_solver_direction_positive(time(i),y(:,i),F, FX);
    k2=fun_17_lorenz_solver_direction_positive(time(i)+step/2,y(:,i)+step*k1/2,F, FX);
    k3=fun_17_lorenz_solver_direction_positive(time(i)+step/2,y(:,i)+step*k2/2,F, FX);
    k4=fun_17_lorenz_solver_direction_positive(time(i)+step,y(:,i)+step*k3,F, FX);
    y(:,i+1)=y(:,i)+step*(k1+2*k2+2*k3+k4)/6;
end
result = y';  % 转为行向量
result(end,[27,28,29])  % 这是结果

% 循环求解
for iter_time = 1:500
    y = y(:, end);
    n=floor((real_b-real_a)/step);       % 步数: 只要步数能够被dot_X覆盖, 就不会有误差
    time(1)=0;             % 时间列表
    % n=1000
    for i=1:n          % 龙格库塔方法进行数值求解
        time(i+1)=time(i)+step;
        k1=fun_17_lorenz_solver_direction_positive(time(i),y(:,i),F, FX);
        k2=fun_17_lorenz_solver_direction_positive(time(i)+step/2,y(:,i)+step*k1/2,F, FX);
        k3=fun_17_lorenz_solver_direction_positive(time(i)+step/2,y(:,i)+step*k2/2,F, FX);
        k4=fun_17_lorenz_solver_direction_positive(time(i)+step,y(:,i)+step*k3,F, FX);
        y(:,i+1)=y(:,i)+step*(k1+2*k2+2*k3+k4)/6;
    end
    result = y';  % 转为行向量
    result(end,[27,28,29])  % 这是结果
end