% 复现Wang 尝试用Yan师兄给的代码进行求解
clear
clc

rng(1)
% 数据准备 - 法1: 无噪声情形, 直接用龙格库塔生成
% step = 0.01;  % 拟合步长为step, 真实轨迹的步长为step/2  (因为:RungeKutta会用到1/2的步长去计算)
% t_u = step;
% b_fit = 100;
% b_test = 20;
% real_b = b_fit+b_test;  % 真实轨迹的结束时间 
% real_a = 0;  % 真实轨迹的起始时间 
% y0 = [-8,7,27];  % 真实轨迹的初始点 
% interp_step = step  % step 或 step/2: #注: 不可以用step/2, 会有很严重的问题
% interp_a = real_a + interp_step;  % 插值F的开始时间  取step/2的中心差分后, 前面会少1个时刻的数据;
% interp_b = real_b - interp_step;  % 插值F的结束时间
% [t,X_n] = fun_8_RungeKutta(@fun_0_lorenz,real_a,interp_step,real_b,y0);  % 求解轨迹: 固定步长: 步长是拟合所用步长的一半
% X = X_n(1:3, 1:10000)';

% 数据准备 - 法2: 有噪声情形
delta=10;b=8/3;r=28;  % Lorenz正常参数
N_fit = 10000
N_test = 2000;
N=N_fit+N_test;
x0 = [-8;7;27];  % 初始值
X_n = [x0];  % R: 容器
t_u = 0.01;  % 单位计算时间
D = 0.01;  % 噪声方差
m = 1;  % R: 哦这个好像是时间窗口之类的
b1=sqrt(2*D*t_u); % 准备使用公式A3
b2=1/2*b1*t_u;
b3=1/sqrt(3)*b2;   
w1=randn(N,3);w2=randn(N,3);  % R: 生成标准白噪声
for i=1:N 
    s1=b1*w1(i,:);  % 
    s2=b2*w1(i,:)+b3*w2(i,:);  % 
    df2=[0;(r-X_n(3,i)-1-X_n(1,i));(X_n(2,i)+X_n(1,i)-b)];  %

    % R: 下面是Lorenz方程, 应该是原方程, 而不是拟合方程.
    fx=delta*(X_n(2,i)-X_n(1,i));  %
    fy=r*X_n(1,i)-X_n(2,i)-X_n(1,i)*X_n(3,i);  %
    fz=X_n(1,i)*X_n(2,i)-b*X_n(3,i);

    % R: 总之应该是往后推一段X_n
    df1=[-delta*fx+delta*fy; (r-X_n(3,i))*fx-fy-X_n(1,i)*fz; X_n(2,i)*fx+X_n(1,i)*fy-b*fz];  % 
    X_t=X_n(:,i)+[fx;fy;fz]*t_u+1/2*t_u*t_u*df1+df2.*s2'+s1';  % R: 用最后一列数据去更新状态, 而且在状态更新时引入了噪声.
    X_n=[X_n,X_t];  % 合并新的列, 最后得到一个3行, 10001列的矩阵. 随着列的增加时间增加.
end
X = X_n(:,1:N_fit)';

% 构建预测矩阵
base = fun_13_baseEuqation_my(X(1:end-1, :));  % (1)构建小base向量, 但不构建最后一行的数据!
BASE = [base*t_u, base*t_u*t_u];  % 这种构造方式并没有完全用复现思路, 仅仅进行了后向预测.
BASE_pinv = pinv(BASE);  %
diff_x = BASE_pinv*(X(2:end,1)-X(1:end-1,1));  % R 做差, 再除以基向量, 这样仿佛能得到, 一个向量使得:这个向量*基向量=与现在的差距.
diff_y = BASE_pinv*(X(2:end,2)-X(1:end-1,2));
diff_z = BASE_pinv*(X(2:end,3)-X(1:end-1,3));

% (失败!)尝试预测: (1)给一个初始点 (2)迭代预测下一个初始点
a_X = X(1,:)
predict_history = [a_X]
for i=1:3000
    a_base = fun_13_baseEuqation_my(a_X);
    a_BASE = [a_base*t_u, a_base*t_u*t_u];
    next_x = a_BASE * diff_x;
    next_y = a_BASE * diff_y;
    next_z = a_BASE * diff_z;
    % 更新数据
    a_X = [next_x, next_y, next_z]
    predict_history = [predict_history; a_X];  % 从这个变量可以看出失败的结果.
end

% 先估算a,b,r
DBASE = [base, 2*t_u*base];  % 原基向量对时间求导
DX = DBASE*diff_x;  % 下一时刻的dot{x}
DY = DBASE*diff_y;  
DZ = DBASE*diff_z;

% 下: 基于lorenz公式, 准备参数之外的项
for i=1:size(DX,1)
    cha_x(i,:)=DX(i,1);  % lorenz 除a的参数项之外的项
end
for i=1:size(DX,1)
    cha_y(i,:)=DY(i,1)+X(i,2)+X(i,1)*X(i,3);  % lorenz 除r的参数项之外的项
end
for i=1:size(DX,1)
    cha_z(i,:)=X(i,1)*X(i,2)-DZ(i,1);  % lorenz 除b的参数项之外的项
end
% 下: 获取各个时刻参数前的系数
for i=1:size(DX,1)  % 
    base_a(i,:)=(X(i,2)-X(i,1));  % lorenz 参数a的系数y-x
end
for i=1:size(DX,1)
    base_b(i,:)=X(i,3);  % lorenz 参数b的系数z
end
for i=1:size(DX,1)
    base_r(i,:)=X(i,1);  % lorenz 参数x的系数r
end
% 下: 准备做除法
base_ap=pinv(base_a);
base_bp=pinv(base_b);
base_rp=pinv(base_r);
% 下: 移项做除法 得到参数回归值
a_value(1)=base_ap*cha_x  % 将系数移至左边, 可得a的值
b_value(1)=base_bp*cha_z  % 同上
r_value(1)=base_rp*cha_y  % 同上


% 尝试预测: 法(1) 基于龙格库塔
% 结果: 在真实数据也是龙格库塔生成的情况下, 大概能往前预测500步左右
% para = [a_value(1), b_value(1), r_value(1)]
% step = 0.01;  % 拟合步长为step, 真实轨迹的步长为step/2  (因为:RungeKutta会用到1/2的步长去计算)
% h = step;
% a = 0;
% b_fit = 0;
% b_test = 20;
% real_b = b_fit + b_test;  % 真实轨迹的结束时间
% real_a = 0;  % 真实轨迹的起始时间
% y0 = X(end,:);  % 真实轨迹的初始点
% n=floor((real_b-real_a)/h);       %步数
% x(1)=a;                 %时间起点
% y=[];
% y(:,1)=y0;              %赋初值，可以是向量，但是要注意维数
% for i=1:n               %龙格库塔方法进行数值求解
%     x(i+1)=x(i)+h;
%     k1=fun_0_lorenz_para(x(i),y(:,i), para);
%     k2=fun_0_lorenz_para(x(i)+h/2,y(:,i)+h*k1/2, para);
%     k3=fun_0_lorenz_para(x(i)+h/2,y(:,i)+h*k2/2, para);
%     k4=fun_0_lorenz_para(x(i)+h,y(:,i)+h*k3, para);
%     y(:,i+1)=y(:,i)+h*(k1+2*k2+2*k3+k4)/6;
% end
% X_pred = y';


% 尝试预测: 法(2) 基于a, b, r的value, 并基于最后一个值往后预测2000步
% 结果: 大概能够往前预测300步左右
delta=a_value(1); b=b_value(1); r=r_value(1);  % Lorenz正常参数
N_test = 2000;
N=N_test;
x0 = X(end,:);  % 初始值
X_pred = [x0'];  % R: 容器
t_u = 0.01;  % 单位计算时间
D = 0;  % 噪声方差
m = 1;  % R: 哦这个好像是时间窗口之类的
b1=sqrt(2*D*t_u); % 准备使用公式A3
b2=1/2*b1*t_u;
b3=1/sqrt(3)*b2;   
w1=randn(N,3);w2=randn(N,3);  % R: 生成标准白噪声
for i=1:N 
    s1=b1*w1(i,:);  % 
    s2=b2*w1(i,:)+b3*w2(i,:);  % 
    df2=[0;(r-X_pred(3,i)-1-X_pred(1,i));(X_pred(2,i)+X_pred(1,i)-b)];  %

    % R: 下面是Lorenz方程, 应该是原方程, 而不是拟合方程.
    fx=delta*(X_pred(2,i)-X_pred(1,i));  %
    fy=r*X_pred(1,i)-X_pred(2,i)-X_pred(1,i)*X_pred(3,i);  %
    fz=X_pred(1,i)*X_pred(2,i)-b*X_pred(3,i);

    % R: 总之应该是往后推一段X_n
    df1=[-delta*fx+delta*fy; (r-X_pred(3,i))*fx-fy-X_pred(1,i)*fz; X_pred(2,i)*fx+X_pred(1,i)*fy-b*fz];  % 
    X_t=X_pred(:,i)+[fx;fy;fz]*t_u+1/2*t_u*t_u*df1+df2.*s2'+s1';  % R: 用最后一列数据去更新状态, 而且在状态更新时引入了噪声.
    X_pred=[X_pred,X_t];  % 合并新的列, 最后得到一个3行, 10001列的矩阵. 随着列的增加时间增加.
end
X_pred = X_pred';

% 画图
figure
plot(X_n(1,:)','linewidth',1);
hold on
t_range = 10000:1:12000;
plot(t_range, X_pred(:,1),'--', 'linewidth',1);
xlim([9500 11000]);


