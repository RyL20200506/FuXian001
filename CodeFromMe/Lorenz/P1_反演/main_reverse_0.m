clear
clc
rng(1);

% 0. 参数设定
real_b = 2000;  % 真实轨迹的结束时间 
real_a = 0;  % 真实轨迹的起始时间 
step = 0.01;  % 拟合步长为step, 真实轨迹的步长为step/2  (因为:RungeKutta会用到1/2的步长去计算)
y0 = [-8,7,27];  % 真实轨迹的初始点 

% 1. 数据准备: 获取已知可观察轨迹
[t,y] = fun_8_RungeKutta(@fun_0_lorenz, real_a, step, real_b, y0); 
real_trajectory = y';
known_X = real_trajectory(:, 1);  % 真实且已知的轨迹
% 验证 - 画图
% figure;
% plot3(real_trajectory(:,1), real_trajectory(:,2), real_trajectory(:,3), 'r')
% grid on;

% 2.1 插值
FX = griddedInterpolant(t, known_X);
% 验证
% known_X(1) - FX(0)
% known_X(10) - FX(0.09)
% known_X(11) - FX(0.1)
% known_X(101) - FX(1)
% known_X(1501) - FX(15)

% 2.2 数据准备: 获取已知可观察轨迹的中心差分
diff_before_X = known_X(2:end-1) - known_X(1:end-2);  % 前向差分: #验证无误
diff_after_X = (known_X(3:end) - known_X(2:end-1));  % 后向差分  #验证无误
diff_center_X = (diff_before_X+diff_after_X)/2;  % 中心差分
dot_X = diff_center_X/(step);  % 数值导数
% 定义插值函数 fun_9_lorenz_solver: 使其能接受外部输入dot_X数据, 并能根据step选取特定的dot_X的数据
time_range = [real_a + step: step: real_b - step];
FdotX = griddedInterpolant(time_range, dot_X);
% 验证 - 原方程做差
% a_t = 100;
% adotx =  10 * ( real_trajectory(a_t,2) -  real_trajectory(a_t,1) );
% FdotX(0.99) - adotx  % 第1个数时间为0, 第2个数时间0.01, 第100个数时间0.99
% a_t = 1000;
% adotx =  10 * ( real_trajectory(a_t,2) -  real_trajectory(a_t,1) );
% FdotX(9.99) - adotx

% 3. 开始求解
n=floor((real_b-real_a)/step);       % 步数: 只要步数能够被dot_X覆盖, 就不会有误差
time(1)=0;             % 时间列表
y = ones(22,1);              % 赋初值，可以是向量，但是要注意维数
% n=1000
for i=1:n          % 龙格库塔方法进行数值求解
    time(i+1)=time(i)+step;
    k1 = fun_lorenz_solver_direction_positive(time(i),y(:,i),FdotX, FX);
    k2 = fun_lorenz_solver_direction_positive(time(i)+step/2,y(:,i)+step*k1/2,FdotX, FX);
    k3 = fun_lorenz_solver_direction_positive(time(i)+step/2,y(:,i)+step*k2/2,FdotX, FX);
    k4 = fun_lorenz_solver_direction_positive(time(i)+step,y(:,i)+step*k3,FdotX, FX);
    y(:,i+1) = y(:,i)+step*(k1+2*k2+2*k3+k4)/6;
end
result = y';  % 转为行向量
result(end,[18,19,20])  % 这是结果

stop;

%% 反向求解: (1)反向沿用时间, 步数
y_n = result(n, :)';              % 沿用之前最后的结果
% n=1000
for i=1:n          % 龙格库塔方法进行数值求解
    index = n - i + 1;
    k1=fun_17_lorenz_solver_direction_negtive(time(index), y_n(:,i), FdotX, FX);
    k2=fun_17_lorenz_solver_direction_negtive(time(index) - step/2, y_n(:,i)+step*k1/2, FdotX, FX);
    k3=fun_17_lorenz_solver_direction_negtive(time(index) - step/2, y_n(:,i)+step*k2/2, FdotX, FX);
    k4=fun_17_lorenz_solver_direction_negtive(time(index) - step, y_n(:,i)+step*k3, FdotX, FX);
    y_n(:,i+1)=y_n(:,i)+step*(k1+2*k2+2*k3+k4)/6;
    y_n(27:29, end)
end
result = y_n';
result(end,[27:29])  % 这是结果

%% 调试
% 首先计算这是原系统在结束时
i = 200000
t = time(i);
Y = result(n, :)';
FXdot = FdotX;

k1 = fun_17_lorenz_solver_direction_positive(time(i),Y, FdotX, FX);  
k1_n = fun_17_lorenz_solver_direction_negtive(time(i),Y, FdotX, FX);  
% 对比上述两个变量会发现, 差异非常大!-1.2078和75.2556; 
% BUT WHY??? 我们关注第22个变量为主, 符号改变, 且扩大了70倍!!!

% 尝试计算k1中的22, 应为 -1.2078
gamma=0.0015;
alpha_m=3;
beta=2;

t = time(i);
dot_x = FdotX(t);
x = FX(t);

% 
hatx = Y(4);
haty = Y(5);
hatz = Y(6);
hatx_hata = Y(7);
hatx_hatb = Y(8);
hatx_hatr = Y(9);
hatx_e1 = Y(10);
hatx_e2 = Y(11);
haty_hata = Y(12);
haty_hatb = Y(13);
haty_hatr = Y(14);
haty_e1 = Y(15);
haty_e2 = Y(16);
hatz_hata = Y(17);
hatz_hatb = Y(18);
hatz_hatr = Y(19);
hatz_e1 = Y(20);
hatz_e2 = Y(21);
D_hata = Y(22);
D_hatb = Y(23);
D_hatr = Y(24);
D_e1 = Y(25);
D_e2 = Y(26);
hata = Y(27);
hatb = Y(28);
hatr = Y(29);
e1 = Y(30);
e2 = Y(31);
% 
direction=1;
delta_a = -alpha_m*D_hata + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ((haty - x) + 2*hata * haty_hata ) 
direction=-1;
delta_a = -alpha_m*D_hata + (-2)*( - dot_x - direction*( hata*(haty-x)) ) * direction * ((haty - x) + 2*hata * haty_hata ) 

% 对于负值


%% 退出调试


% 画出参数的运动图
figure
plot(result(:,25), 'LineWidth',1.5) % a
hold on
plot(result(:,26), 'LineWidth',1.5) % b
hold on 
plot(result(:,27), 'LineWidth',1.5) % r
hold on
legend('a','b','r') ;
xlabel('\it Step \rm', 'fontsize',17);
ylabel('\it Value \rm', 'fontsize',17);
% 自定义刻度
xtickformat('%.1f');
ytickformat('%.1f');
ax = gca;
ax.XAxis.Exponent = 5;
ylim([0 30])
% 调整字体大小
set(gca,'FontSize',14)  %是设置刻度字体大小
% % 画出参数的运动图

% figure
% plot(result(:,25)) % a
% hold on 
% plot(result(:,26)) % b
% hold on 
% plot(result(:,27)) % r
% hold on 
% legend('a','b','r') 

stop  % 下面都是调试

% 画图(1): 不要遮挡有效数据
% 画图(2): 放在左边
