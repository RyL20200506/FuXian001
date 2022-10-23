% ����Yan. ʹ���Լ�ʵ�ֵ��������, ȥ��ϲ�ַ�, �������
clear
clc

step = 0.01;  % ��ϲ���Ϊstep, ��ʵ�켣�Ĳ���Ϊstep/2  (��Ϊ:RungeKutta���õ�1/2�Ĳ���ȥ����)
real_b = 2000;  % ��ʵ�켣�Ľ���ʱ�� 
real_a = 0;  % ��ʵ�켣����ʼʱ�� 
y0 = [-8,7,27];  % ��ʵ�켣�ĳ�ʼ�� 
interp_step = step  % step �� step/2: #ע: ��������step/2, ���к����ص�����
interp_a = real_a + interp_step;  % ��ֵF�Ŀ�ʼʱ��  ȡstep/2�����Ĳ�ֺ�, ǰ�����1��ʱ�̵�����;
interp_b = real_b - interp_step;  % ��ֵF�Ľ���ʱ��

% ����׼��: ��ȡ��֪�ɹ۲�켣
[t,real_trajectory] = fun_8_RungeKutta(@fun_0_lorenz,real_a,interp_step,real_b,y0);  % ���켣: �̶�����: ������������ò�����һ��
real_trajectory = real_trajectory';
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
y=[-8; 7; 27; ones(26,1)];              % ����ֵ������������������Ҫע��ά��

% n=1000
for i=1:n          % �����������������ֵ���
    time(i+1)=time(i)+step;
    k1=fun_9_lorenz_solver(time(i),y(:,i),F, FX);
    k2=fun_9_lorenz_solver(time(i)+step/2,y(:,i)+step*k1/2,F, FX);
    k3=fun_9_lorenz_solver(time(i)+step/2,y(:,i)+step*k2/2,F, FX);
    k4=fun_9_lorenz_solver(time(i)+step,y(:,i)+step*k3,F, FX);
    y(:,i+1)=y(:,i)+step*(k1+2*k2+2*k3+k4)/6;
end
result = y';
result(end,[25,26,27])  % ���ǽ��


% �����������˶�ͼ
figure
plot(result(:,25)) % a
hold on 
plot(result(:,26)) % b
hold on 
plot(result(:,27)) % r
hold on 
legend('a','b','r') 

stop  % ���涼�ǵ���

% ������ϲ���Ĺ켣ͼ
figure;
plot3(result(:,1),result(:,2),result(:,3),'r')
grid on;

% ����ԭʼ����ʵ�켣ͼ
figure
plot3(real_trajectory(:,1),real_trajectory(:,2),real_trajectory(:,3));  % ��ͼ
grid on;  % ������

% ����켣���
save('mat_9_result','result')

% Ŀǰ������: F��ֵ(3)��΢�ַ��̵���ʵ������һ��, ��F��ֵ��dot_X�Ƕ�Ӧ��.
% (1)����: �Ա����ֲ�ͬ��ʵ�켣, 3������, ��29�������Ĺ켣�����Ƿ�һ��-> �������������ǵĽ����ȫһ��
[t,real_3] = fun_8_RungeKutta(@fun_9_lorenz,real_a,step/2,real_b,[-8,7,27]);  % ���켣: 3��Ԫ�صĹ켣
real_3 = real_3';
known_X_3 = real_3(:, 1);  % ��ʵ����֪�Ĺ켣
[t,real_29] = fun_8_RungeKutta(@fun_9_lorenz_solver,real_a,step/2,real_b,[-8; 7; 27; ones(26,1)]);  % ���켣: 29��Ԫ�صĹ켣
real_29 = real_29';
known_X_29 = real_29(:, 1);  % ��������ȫһ��.

% ����: �ԱȲ�ַ��õ��ĵ���, ����ʵ�켣�ĵ����Ƿ�һ��: �о����ܲ�һ��
fun_8_RungeKutta(@fun_9_lorenz,real_a,step/2,real_b,y0);

% (2)����: �ڲ�������ͬ�������, ��dotX�����Ӱ��; ��Ϊ��ʵ���ݵĲ�����step/2, �������ݵĲ�����step 
% -> ȷʵ��������⵼�µ�! �������dotX��ɺܴ��Ӱ��, ���ǵĹ켣Ӧ�ö��кܴ�����!
% -> Q:��ô, dot_XӦ����εõ���?
[t,real_3] = fun_8_RungeKutta(@fun_9_lorenz,real_a, step/2, real_b,[-8,7,27]);  % ���켣: 3��Ԫ�صĹ켣
real_3 = real_3';
known_X_3 = real_3(:, 1);  % ��ʵ����֪�Ĺ켣
diff_center_X_3 = fun_9_diff_center(known_X_3);
dot_X_3 = diff_center_X_3/(step/2);
time_range = [real_a+step/2: step/2: real_b-step/2];
F_3 = griddedInterpolant(time_range, dot_X_3);

[t,real_29] = fun_8_RungeKutta(@fun_9_lorenz_solver,real_a, step, real_b,[-8; 7; 27; ones(26,1)]);  % ���켣: 29��Ԫ�صĹ켣
real_29 = real_29';
known_X_29 = real_29(:, 1);  % ��ʵ����֪�Ĺ켣
diff_center_X_29 = fun_9_diff_center(known_X_29);
dot_X_29 = diff_center_X_29/(step);
time_range = [real_a+step: step: real_b-step];
F_29 = griddedInterpolant(time_range, dot_X_29);

% (3)����: �ɴ˿�֪, ���ǽ�����ͬ����, �������Բ�ֵ�ķ���ȥ���
% -> 