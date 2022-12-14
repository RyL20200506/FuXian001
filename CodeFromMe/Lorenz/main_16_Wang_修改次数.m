% 用于验证:
% Xie师姐的代码效果好的原因, 是由于baseEquation1展开到了5维所导致的
% 所以:(1)将baseEquation修改为3维, (2)将baseEquation修改为10维

% 效果: 对于不同方法生成的数据, Xie师姐的代码都有很好的兼容性!

close all;clear;clc;

% 生成待预测数据
[t1,h1]=fun_8_RungeKutta(@fun_0_lorenz, 0, 0.01, 200, [12 4 0]);
% plot3(h1(1,:),h1(2,:),h1(3,:),'r'); grid on;
history = h1';
save('mat_14_history._from_me.mat', 'history')

% 准备数据
history_train = history(1:10000, :);
predict_point = history_train(end, :);
predict_length = 5000;

% 开始预测
% future_values = fun_15_GeneralPredict_Xie(history_train, predict_point, predict_length);
future_values = fun_16_GeneralPredict_by_para_10(history_train, predict_point, predict_length);

% 画预测的图: Xie师姐的方法很好!
figure
history_test = history(10000:end, :)
plot(history_test(:,1))
hold on
plot(future_values(:,1))
xlim([0,3000])

stop;

%% 更换数据集继续
close all;clear;clc;

% 生成待预测数据
load mat_15_history_from_yan.mat;

% 准备数据
history_train = history(1:10000, :);
predict_point = history_train(end, :);
predict_length = 5000;

% 开始预测
future_values = fun_15_GeneralPredict_Xie(history_train, predict_point, predict_length);

% 画预测的图: Xie师姐的方法很好!
figure
history_test = history(10000:end, :)
plot(history_test(:,1))
hold on
plot(future_values(:,1))

stop;

