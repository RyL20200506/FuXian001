
% ���ڻ���base�������õ�����ֱ���õ�B����
% ����:�������Ӽ���
function [B] = fun_13_1_BuildB(base_all, t_unit)

    % ��������
%     base_all = base;  % �������л����ɵ��˱���
%     t_unit = 0.01;

    [num_row,num_column] = size(base_all);  % ��ȡ�������Ĺ��

    % ������ѭ��
    B = zeros(3*(num_row-2), num_column*2);  % ��Ϊ��t��t^2����ά��, ����Ҫb*2
%     B = [];
    for i=1:num_row
        a_base_row = base_all(i,:);  % ȡ����ǰ��base����
        % ���ɵ�ǰbase��������Ӧ��3��
        row_last = [a_base_row*(-t_unit), a_base_row*(t_unit)^2];
        row_present = zeros(1, size(row_last,2));
        row_next = [a_base_row*t_unit, a_base_row*(t_unit)^2];
        % ��3����ӽ�ȥ
        B(i:i+2, :) = [row_last; row_present; row_next];
    end
end


