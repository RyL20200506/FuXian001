% ��Ҫ������main12, ����main1���, ���ǿ�������һ����main1��ģ��(��Ϊ�⻹�����ȷ�����)
% �����ܲ����ñȽ��ٵĲ�������

% ����Yan. ʹ���Լ�ʵ�ֵ��������, ȥ��ϲ�ַ�, �������
clear
clc
rng(1);
step = 0.01;  % ��ϲ���Ϊstep, ��ʵ�켣�Ĳ���Ϊstep/2  (��Ϊ:RungeKutta���õ�1/2�Ĳ���ȥ����)
real_b = 2000;  % ��ʵ�켣�Ľ���ʱ�� 
real_a = 0;  % ��ʵ�켣����ʼʱ�� 
y0 = [-8,7,27];  % ��ʵ�켣�ĳ�ʼ�� 
interp_step = step  % step �� step/2: #ע: ��������step/2, ���к����ص�����
interp_a = real_a + interp_step;  % ��ֵF�Ŀ�ʼʱ��  ȡstep/2�����Ĳ�ֺ�, ǰ�����1��ʱ�̵�����;
interp_b = real_b - interp_step;  % ��ֵF�Ľ���ʱ��

% ����׼��: ��ȡ��֪�ɹ۲�켣
[t,y] = fun_8_RungeKutta(@fun_11_lorenz_noise, real_a, interp_step, real_b, y0); 
real_trajectory = y';
known_X = real_trajectory(:, 1);  % ��ʵ����֪�Ĺ켣
% ��ֵ: ֱ�Ӳ�X��ֵ
FX = griddedInterpolant(t, known_X);

% ����׼��: ��ȡ��֪�ɹ۲�켣�����Ĳ��, �Ӷ�
diff_before_X = known_X(2:end-1) - known_X(1:end-2);  % ǰ����: #��֤����
diff_after_X = (known_X(3:end) - known_X(2:end-1));  % ������  #��֤����
diff_center_X = (diff_before_X+diff_after_X)/2;  % ���Ĳ��
dot_X = diff_center_X/(interp_step);  % ��ֵ����
save('mat_9_dot_X', 'dot_X');  % �����м����

% �����ֵ���� fun_9_lorenz_solver: ʹ���ܽ����ⲿ����dot_X����, ���ܸ���stepѡȡ�ض���dot_X������
time_range = [interp_a: interp_step: interp_b];
F = griddedInterpolant(time_range, dot_X);

% ��ʼ���
n=floor((real_b-real_a)/step);       % ����: ֻҪ�����ܹ���dot_X����, �Ͳ��������
time(1)=0;             % ʱ���б�
% y = [-8; 7; 27; ones(28,1)];              % ����ֵ������������������Ҫע��ά��
y = [-8; 7; 27; ones(26,1)];              % ����ֵ������������������Ҫע��ά��
% n=1000
for i=1:n          % �����������������ֵ���
    time(i+1)=time(i)+step;
    k1=fun_17_lorenz_solver_direction_positive_o(time(i),y(:,i),F, FX);
    k2=fun_17_lorenz_solver_direction_positive_o(time(i)+step/2,y(:,i)+step*k1/2,F, FX);
    k3=fun_17_lorenz_solver_direction_positive_o(time(i)+step/2,y(:,i)+step*k2/2,F, FX);
    k4=fun_17_lorenz_solver_direction_positive_o(time(i)+step,y(:,i)+step*k3,F, FX);
    y(:,i+1)=y(:,i)+step*(k1+2*k2+2*k3+k4)/6;
end
result = y';  % תΪ������
result(end,[25,26,27])  % ���ǽ��
% result(end,[27,28,29])  % ���ǽ��

stop;

% �������: (1)��������ʱ��, ����
y_n = result(end, :)';              % ����֮ǰ���Ľ��
% n=1000
for i=1:n          % �����������������ֵ���
    index = n - i + 2;
    k1=fun_17_lorenz_solver_direction_negtive(time(index), y_n(:,i), F, FX);
    k2=fun_17_lorenz_solver_direction_negtive(time(index) - step/2, y_n(:,i)+step*k1/2, F, FX);
    k3=fun_17_lorenz_solver_direction_negtive(time(index) - step/2, y_n(:,i)+step*k2/2, F, FX);
    k4=fun_17_lorenz_solver_direction_negtive(time(index) - step, y_n(:,i)+step*k3, F, FX);
    y_n(:,i+1)=y_n(:,i)+step*(k1+2*k2+2*k3+k4)/6;
    y_n(25:27 ,end)
end
result = y_n';
% result(end,[25,26,27])  % ���ǽ��
result(end,[25,26,27])  % ���ǽ��


% �����������˶�ͼ
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
% �Զ���̶�
xtickformat('%.1f');
ytickformat('%.1f');
ax = gca;
ax.XAxis.Exponent = 5;
ylim([0 30])
% ���������С
set(gca,'FontSize',14)  %�����ÿ̶������С
% % �����������˶�ͼ

% figure
% plot(result(:,25)) % a
% hold on 
% plot(result(:,26)) % b
% hold on 
% plot(result(:,27)) % r
% hold on 
% legend('a','b','r') 

stop  % ���涼�ǵ���

% ��ͼ(1): ��Ҫ�ڵ���Ч����
% ��ͼ(2): �������
 