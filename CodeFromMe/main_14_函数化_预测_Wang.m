%% HOPE
% 在(1)不知道函数形式, (2)只输入一堆数据的情况下, (3)做出比较准确的预测
close all;clear;clc;

% 生成待预测数据
[t1,h1]=fun_8_RungeKutta(@fun_0_lorenz, 0, 0.01, 200, [12 4 0]);
plot3(h1(1,:),h1(2,:),h1(3,:),'r'); grid on;
history = h1';
save('mat_14_history.mat', 'history')

% 准备数据
history_train = history(1:10000, :);
predict_point = history_train(end, :);
predict_length = 5000;

% 开始预测
future_values = fun_14_GeneralPredict_Wang(history_train, predict_point, predict_length);

% 画预测的图: 效果一般, 大概预测
history_test = history(10000:end, :)
plot(history_test(:,1))
hold on
plot(future_values(:,1))

stop;

%% 更换数据: 附录A3生成数据
close all;clear;clc;
delta=10;b=8/3;r=28;
N=15000;
% x0=[-1;3;4];
x0=[-8;7;27];
X_n=[x0];
t0=0.01;%单位计算时间
m=1;
D=0;%噪声方差
% D=0;
b1=sqrt(2*D*t0);b2=1/2*b1*t0;b3=1/sqrt(3)*b2;
w1=randn(N,3);w2=randn(N,3);
for i=1:N
    s1=b1*w1(i,:);
    s2=b2*w1(i,:)+b3*w2(i,:);
    df2=[0;(r-X_n(3,i)-1-X_n(1,i));(X_n(2,i)+X_n(1,i)-b)];
    fx=delta*(X_n(2,i)-X_n(1,i));
    fy=r*X_n(1,i)-X_n(2,i)-X_n(1,i)*X_n(3,i);
    fz=X_n(1,i)*X_n(2,i)-b*X_n(3,i);
    df1=[-delta*fx+delta*fy;(r-X_n(3,i))*fx-fy-X_n(1,i)*fz;X_n(2,i)*fx+X_n(1,i)*fy-b*fz]; 
    X_t=X_n(:,i)+[fx;fy;fz]*t0+1/2*t0*t0*df1+df2.*s2'+s1';
    X_n=[X_n,X_t];
end
history = X_n';

% 准备数据
history_train = history(1:10000, :);
predict_point = history_train(end, :);
predict_length = 4000;

% 开始预测
future_values = fun_14_GeneralPredict_Wang(history_train, predict_point, predict_length);

% 画预测的图: 预测范围2800
figure
history_test = history(10000:end, :)
plot(history_test(:,1))
hold on
plot(future_values(:,1))

stop;

%% 更换为xie的数据集
close all;clear;clc;

% 生成待预测数据
load mat_14_history_from_xie.mat;

% 准备数据
history_train = history(1:8000, :);
start_from = 1;
predict_point = history_train(start_from, :);
predict_length = 2000;

% 开始预测
future_values = fun_14_GeneralPredict_Wang(history_train, predict_point, predict_length);

% 画预测的图
figure
history_test = history(start_from:end, :)
plot(history_test(:,1))
hold on
plot(future_values(:,1))

stop;








