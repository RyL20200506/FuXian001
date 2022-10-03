% ����Yan. ��main_10�����Ͻ������޸�Ϊ, ÿһ������������Ӱ�쵽��һʱ�̵Ľ��ֵ.
clear
clc

step = 0.01;  % ��ϲ���Ϊstep, ��ʵ�켣�Ĳ���Ϊstep/2  (��Ϊ:RungeKutta���õ�1/2�Ĳ���ȥ����)
real_b = 20;  % ��ʵ�켣�Ľ���ʱ�� 
real_a = 0;  % ��ʵ�켣����ʼʱ�� 
y0 = [-8,7,27];  % ��ʵ�켣�ĳ�ʼ��
interp_step = step  % step �� step/2: #ע: ��������step/2, ���к����ص�����
interp_a = real_a + interp_step;  % ��ֵF�Ŀ�ʼʱ��  ȡstep/2�����Ĳ�ֺ�, ǰ�����1��ʱ�̵�����;
interp_b = real_b - interp_step;  % ��ֵF�Ľ���ʱ��
D = 0.2;

% ����׼��: ��ȡ��֪�ɹ۲�켣, �������������汾
a=real_a;
b=real_b;
h=interp_step;
n=floor((b-a)/h);       %����
x(1)=a;                 %ʱ�����
y(:,1)=y0;              %����ֵ������������������Ҫע��ά��
rng(1);                 % ���������������
for i=1:n               %�����������������ֵ���
    x(i+1)=x(i)+h;
    k1=fun_0_lorenz(x(i),y(:,i));
    k2=fun_0_lorenz(x(i)+h/2,y(:,i)+h*k1/2);
    k3=fun_0_lorenz(x(i)+h/2,y(:,i)+h*k2/2);
    k4=fun_0_lorenz(x(i)+h,y(:,i)+h*k3);
    y(:,i+1)=y(:,i)+h*(k1+2*k2+2*k3+k4)/6 + D*randn(1);
end
real_trajectory = y';

% ����ԭʼ����ʵ�켣ͼ
figure
plot3(real_trajectory(:,1),real_trajectory(:,2),real_trajectory(:,3));  % ��ͼ
grid on;  % ������

% ֻ��x-z���ͼ��
plot(real_trajectory(:,1),real_trajectory(:,3));  % xz��ͼ�켣

% ע: ���̵�˵, ���滭��ͼ��, ������������ϣ��������ͼ��; ��Ҫԭ������, ԭ������dot{x}�ϼ�����, ���������ڽ���ϼ�����
% ��Ҫ�������������дһ����������Lorenz����

% �ù�ʽ(1)ȥ������.
[t,y] = fun_8_RungeKutta(@fun_11_lorenz_noise,1,interp_step,100,y0); 
real_trajectory = y';
figure
plot(real_trajectory(:,1),real_trajectory(:,3));  % xz��ͼ�켣
