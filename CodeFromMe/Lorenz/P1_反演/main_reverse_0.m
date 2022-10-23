clear
clc
rng(1);

% 0. �����趨
real_b = 2000;  % ��ʵ�켣�Ľ���ʱ�� 
real_a = 0;  % ��ʵ�켣����ʼʱ�� 
step = 0.01;  % ��ϲ���Ϊstep, ��ʵ�켣�Ĳ���Ϊstep/2  (��Ϊ:RungeKutta���õ�1/2�Ĳ���ȥ����)
y0 = [-8,7,27];  % ��ʵ�켣�ĳ�ʼ�� 

% 1. ����׼��: ��ȡ��֪�ɹ۲�켣
[t,y] = fun_8_RungeKutta(@fun_0_lorenz, real_a, step, real_b, y0); 
real_trajectory = y';
known_X = real_trajectory(:, 1);  % ��ʵ����֪�Ĺ켣
% ��֤ - ��ͼ
% figure;
% plot3(real_trajectory(:,1), real_trajectory(:,2), real_trajectory(:,3), 'r')
% grid on;

% 2.1 ��ֵ
FX = griddedInterpolant(t, known_X);
% ��֤
% known_X(1) - FX(0)
% known_X(10) - FX(0.09)
% known_X(11) - FX(0.1)
% known_X(101) - FX(1)
% known_X(1501) - FX(15)

% 2.2 ����׼��: ��ȡ��֪�ɹ۲�켣�����Ĳ��
diff_before_X = known_X(2:end-1) - known_X(1:end-2);  % ǰ����: #��֤����
diff_after_X = (known_X(3:end) - known_X(2:end-1));  % ������  #��֤����
diff_center_X = (diff_before_X+diff_after_X)/2;  % ���Ĳ��
dot_X = diff_center_X/(step);  % ��ֵ����
% �����ֵ���� fun_9_lorenz_solver: ʹ���ܽ����ⲿ����dot_X����, ���ܸ���stepѡȡ�ض���dot_X������
time_range = [real_a + step: step: real_b - step];
FdotX = griddedInterpolant(time_range, dot_X);
% ��֤ - ԭ��������
% a_t = 100;
% adotx =  10 * ( real_trajectory(a_t,2) -  real_trajectory(a_t,1) );
% FdotX(0.99) - adotx  % ��1����ʱ��Ϊ0, ��2����ʱ��0.01, ��100����ʱ��0.99
% a_t = 1000;
% adotx =  10 * ( real_trajectory(a_t,2) -  real_trajectory(a_t,1) );
% FdotX(9.99) - adotx

% 3. ��ʼ���
n=floor((real_b-real_a)/step);       % ����: ֻҪ�����ܹ���dot_X����, �Ͳ��������
time(1)=0;             % ʱ���б�
y = ones(22,1);              % ����ֵ������������������Ҫע��ά��
% n=1000
for i=1:n          % �����������������ֵ���
    time(i+1)=time(i)+step;
    k1 = fun_lorenz_solver_direction_positive(time(i),y(:,i),FdotX, FX);
    k2 = fun_lorenz_solver_direction_positive(time(i)+step/2,y(:,i)+step*k1/2,FdotX, FX);
    k3 = fun_lorenz_solver_direction_positive(time(i)+step/2,y(:,i)+step*k2/2,FdotX, FX);
    k4 = fun_lorenz_solver_direction_positive(time(i)+step,y(:,i)+step*k3,FdotX, FX);
    y(:,i+1) = y(:,i)+step*(k1+2*k2+2*k3+k4)/6;
end
result = y';  % תΪ������
result(end,[18,19,20])  % ���ǽ��

stop;

%% �������: (1)��������ʱ��, ����
y_n = result(n, :)';              % ����֮ǰ���Ľ��
% n=1000
for i=1:n          % �����������������ֵ���
    index = n - i + 1;
    k1=fun_17_lorenz_solver_direction_negtive(time(index), y_n(:,i), FdotX, FX);
    k2=fun_17_lorenz_solver_direction_negtive(time(index) - step/2, y_n(:,i)+step*k1/2, FdotX, FX);
    k3=fun_17_lorenz_solver_direction_negtive(time(index) - step/2, y_n(:,i)+step*k2/2, FdotX, FX);
    k4=fun_17_lorenz_solver_direction_negtive(time(index) - step, y_n(:,i)+step*k3, FdotX, FX);
    y_n(:,i+1)=y_n(:,i)+step*(k1+2*k2+2*k3+k4)/6;
    y_n(27:29, end)
end
result = y_n';
result(end,[27:29])  % ���ǽ��

%% ����
% ���ȼ�������ԭϵͳ�ڽ���ʱ
i = 200000
t = time(i);
Y = result(n, :)';
FXdot = FdotX;

k1 = fun_17_lorenz_solver_direction_positive(time(i),Y, FdotX, FX);  
k1_n = fun_17_lorenz_solver_direction_negtive(time(i),Y, FdotX, FX);  
% �Ա��������������ᷢ��, ����ǳ���!-1.2078��75.2556; 
% BUT WHY??? ���ǹ�ע��22������Ϊ��, ���Ÿı�, ��������70��!!!

% ���Լ���k1�е�22, ӦΪ -1.2078
gamma=0.0015;
alpha_m=3;
beta=2;

t = time(i);
dot_x = FdotX(t);
x = FX(t);

% 
hatx = Y(4);
haty = Y(5);
hatz = Y(6);
hatx_hata = Y(7);
hatx_hatb = Y(8);
hatx_hatr = Y(9);
hatx_e1 = Y(10);
hatx_e2 = Y(11);
haty_hata = Y(12);
haty_hatb = Y(13);
haty_hatr = Y(14);
haty_e1 = Y(15);
haty_e2 = Y(16);
hatz_hata = Y(17);
hatz_hatb = Y(18);
hatz_hatr = Y(19);
hatz_e1 = Y(20);
hatz_e2 = Y(21);
D_hata = Y(22);
D_hatb = Y(23);
D_hatr = Y(24);
D_e1 = Y(25);
D_e2 = Y(26);
hata = Y(27);
hatb = Y(28);
hatr = Y(29);
e1 = Y(30);
e2 = Y(31);
% 
direction=1;
delta_a = -alpha_m*D_hata + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ((haty - x) + 2*hata * haty_hata ) 
direction=-1;
delta_a = -alpha_m*D_hata + (-2)*( - dot_x - direction*( hata*(haty-x)) ) * direction * ((haty - x) + 2*hata * haty_hata ) 

% ���ڸ�ֵ


%% �˳�����


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
