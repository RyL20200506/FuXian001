% 希望实现: 
% 在(1)不知道函数形式, (2)只输入一堆数据的情况下, (3)做出准确预测
function [future_values] = fun_15_GeneralPredict_Xie(history, predict_start_point, predict_length)
    % 输入参数: history - 历史数据
    % 输入参数: predict_start_point - 预测起始点
    % 输入参数: predict_length - 预测长度
    
    % 测试输入
%     close all;clear;clc;
%     load mat_14_history.mat;  % 得到一个history的变量
%     history = history(1:10000, :);
%     predict_start_point = history(end, :);
%     predict_length = 5000;
    t0 = 0.01;
    
    % 循环构建BASE向量
%     [history_length, ~] = size(history); 
%     BASE = [];  % 所有时刻的基向量
%     A_VAL = [];  % 每个时刻
%     m = 1;  % 时间窗口
%     for i = m+1:history_length-m
%         for j = -m:m
%             a_history = history(i, :);
%             a_base = fun_15_baseEuqation1(a_history);
%             BASE = [BASE; a_base * j * t0, a_base * j^2 * t0^2];  % 希望通过这里的BASE
%             A_VAL = [A_VAL; history(i+j, :) - history(i, :)];  % 来预测这里的A_VAL
%         end
%     end
    
    % 用其他的方式获取BASE和A_VAL
    A_VAL = history(2:end, :)-history(1:end-1,:);
    base = fun_15_baseEuqation1(history(1:end-1,:))';
    BASE = [base*t0, base*t0^2];
    
    % 参数拟合
    BASE_inv = pinv(BASE);
    diff = BASE_inv * A_VAL;
    
    % 开始预测
    future_values = [predict_start_point];
    for i = 1: predict_length
        a_future = future_values(i, :);
        a_base = fun_15_baseEuqation1(a_future)';
        a_BASE = [a_base * t0, a_base * t0^2];
        next_values = a_future + a_BASE * diff;
        future_values = [future_values; next_values];
    end
end


