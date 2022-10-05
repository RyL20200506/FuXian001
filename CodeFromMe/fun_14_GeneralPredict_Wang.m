% ϣ��ʵ��: 
% ��(1)��֪��������ʽ, (2)ֻ����һ�����ݵ������, (3)����׼ȷԤ��
function [future_values] = fun_14_GeneralPredict_Wang(history, predict_start_point, predict_length)
    % �������: history - ��ʷ����
    % �������: predict_start_point - Ԥ����ʼ��
    % �������: predict_length - Ԥ�ⳤ��
    
    % ��������
%     load mat_14_history.mat;  % �õ�һ��history�ı���
%     history = history(1:10000, :);
%     predict_start_point = history(end, :);
%     predict_length = 5000;
    
    % ѭ������BASE����
    [history_length, ~] = size(history); 
    BASE = [];  % ����ʱ�̵Ļ�����
    A_VAL = [];  % ÿ��ʱ��
    m = 1;  % ʱ�䴰��
    t0 = 0.01;
    for i = m+1:history_length-m
        for j = -m:m
            a_history = history(i, :);
            a_base = fun_15_baseEuqation1(a_history);
            BASE = [BASE; a_base * j * t0, a_base * j^2 * t0^2];  % ϣ��ͨ�������BASE
            A_VAL = [A_VAL; history(i+j, :) - history(i, :)];  % ��Ԥ�������A_VAL
        end
    end
    
    % �������
    BASE_pinv = pinv(BASE);
    diff = BASE_pinv * A_VAL;
    
    % ��ʼԤ��
    future_values = [predict_start_point];
    for i = 1: predict_length
        a_future = future_values(i, :);
        a_base = fun_13_baseEuqation(a_future);
        a_BASE = [a_base * j * t0, a_base * j^2 * t0^2];
        next_values = a_future + a_BASE * diff;
        future_values = [future_values; next_values];
    end
end

