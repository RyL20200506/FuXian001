% 复现: Wang - 重构动力学
% 函数: 输入参数, 初值, 目标时间, 根据多项式来重构整条轨迹

function data_history = fun_3_4_predict(PARA, X0, T)

%测试输入
% PARA = ones(48,1);
% X0 = [1,1,1];
% T=1;  % 默认仿真次数

t_unit = 0.01;  % 原文的预测时长: 所以真实时间=T*t_unit
data_history = X0;  % 预测时间序列容器, 行数随着时间增加, 而列为不同维度

for time=1:T
    % 获取最新一行的状态变量
    x = data_history(end,1);
    y = data_history(end,2);
    z = data_history(end,3);
    
    % 获取参数
    a_x=PARA(1);
    b_x_x=PARA(2);
    b_x_y=PARA(3);
    b_x_z=PARA(4);
    c_x_x_y=PARA(5);
    c_x_x_z=PARA(6);
    c_x_y_z=PARA(7);
    d_x_x_y_z=PARA(8);
    a1_x=PARA(9);
    b1_x_x=PARA(10);
    b1_x_y=PARA(11);
    b1_x_z=PARA(12);
    c1_x_x_y=PARA(13);
    c1_x_x_z=PARA(14);
    c1_x_y_z=PARA(15);
    d1_x_x_y_z=PARA(16);
    a_y=PARA(17);
    b_y_x=PARA(18);
    b_y_y=PARA(19);
    b_y_z=PARA(20);
    c_y_x_y=PARA(21);
    c_y_x_z=PARA(22);
    c_y_y_z=PARA(23);
    d_y_x_y_z=PARA(24);
    a1_y=PARA(25);
    b1_y_x=PARA(26);
    b1_y_y=PARA(27);
    b1_y_z=PARA(28);
    c1_y_x_y=PARA(29);
    c1_y_x_z=PARA(30);
    c1_y_y_z=PARA(31);
    d1_y_x_y_z=PARA(32);
    a_z=PARA(33);
    b_z_x=PARA(34);
    b_z_y=PARA(35);
    b_z_z=PARA(36);
    c_z_x_y=PARA(37);
    c_z_x_z=PARA(38);
    c_z_y_z=PARA(39);
    d_z_x_y_z=PARA(40);
    a1_z=PARA(41);
    b1_z_x=PARA(42);
    b1_z_y=PARA(43);
    b1_z_z=PARA(44);
    c1_z_x_y=PARA(45);
    c1_z_x_z=PARA(46);
    c1_z_y_z=PARA(47);
    d1_z_x_y_z=PARA(48);

    % 计算: 中间变量
    g_x = a_x + b_x_x*x+b_x_y*y+b_x_z*z ...
        + c_x_x_y*x*y+c_x_x_z*x*z+c_x_y_z*y*z ...
        + d_x_x_y_z*x*y*z;
    h_x = a1_x + b1_x_x*x+b1_x_y*y+b1_x_z*z ...
        + c1_x_x_y*x*y+c1_x_x_z*x*z+c1_x_y_z*y*z...
        + d1_x_x_y_z*x*y*z;
    g_y = a_y + b_y_x*x+b_y_y*y+b_y_z*z...
        + c_y_x_y*x*y+c_y_x_z*x*z+c_y_y_z*y*z...
        + d_y_x_y_z*x*y*z;
    h_y = a1_y + b1_y_x*x+b1_y_y*y+b1_y_z*z...
        + c1_y_x_y*x*y+c1_y_x_z*x*z+c1_y_y_z*y*z...
        + d1_y_x_y_z*x*y*z;
    g_z = a_z + b_z_x*x+b_z_y*y+b_z_z*z...
        + c_z_x_y*x*y+c_z_x_z*x*z+c_z_y_z*y*z...
        + d_z_x_y_z*x*y*z;
    h_z = a1_z + b1_z_x*x+b1_z_y*y+b1_z_z*z...
        + c1_z_x_y*x*y+c1_z_x_z*x*z+c1_z_y_z*y*z...
        + d1_z_x_y_z*x*y*z;

    % 计算: 下一刻的数值
    next_x = x + g_x*t_unit + h_x*t_unit*t_unit;
    next_y = y + g_y*t_unit + h_y*t_unit*t_unit;
    next_z = z + g_z*t_unit + h_z*t_unit*t_unit;

    % 记录
    data_history = [data_history; next_x, next_y, next_z];
end
data_history = data_history(2:end,:);  % 这里我们会去掉刚开始的首行数据, 只输出预测数据

% 该预测函数会出现的一个问题是, 当T过长, 轨迹可能会发散出去. 这不是我们希望的.
% 我们可以构造出更加精确的数据集, 就是只看下一时刻的预测, 然后在当前状态已知的情况下, 希望下一时刻越接近越好.








