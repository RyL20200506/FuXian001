
% 用于基于base基向量得到可以直接用的B向量
% 方法:逐行增加即可
function [B] = fun_13_1_BuildB(base_all, t_unit)

    % 测试输入
%     base_all = base;  % 主程序中会生成的人变量
%     t_unit = 0.01;

    [num_row,num_column] = size(base_all);  % 获取基向量的规格

    % 对行数循环
    B = zeros(3*(num_row-2), num_column*2);  % 因为有t和t^2两个维度, 所以要b*2
%     B = [];
    for i=1:num_row
        a_base_row = base_all(i,:);  % 取出当前的base向量
        % 生成当前base向量所对应的3行
        row_last = [a_base_row*(-t_unit), a_base_row*(t_unit)^2];
        row_present = zeros(1, size(row_last,2));
        row_next = [a_base_row*t_unit, a_base_row*(t_unit)^2];
        % 将3行添加进去
        B(i:i+2, :) = [row_last; row_present; row_next];
    end
end


