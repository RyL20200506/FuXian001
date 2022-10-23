% ���Ի�����lorenz��ͼ��
close all
clear
clc


%% ��ԭʧ��
y0 = [1,1,1]; 
% �������: ֱ�ӵ���
length = 2;
[t,y] = fun_8_RungeKutta(@fun_0_lorenz,0,0.01,length,y0);
y = y';

% ��ͼ
figure
plot(y(:,1),y(:,3));  % xz��ͼ�켣
grid on;


% ���� ������� ���ú���
% ���һ����
the_point = y(end, :);
[t_n, y_n] = fun_8_RungeKutta(@fun_18_negative_lorenz, 0, 0.01, length, the_point);
y_n = y_n';

figure
back_n = length/0.01
plot(y_n(1:back_n ,1), y_n(1:back_n ,3));  % xz��ͼ�켣
grid on;

stop
%% �ɹ���ԭ
close all
clear
clc

y0 = [1,1,1]; 
length = 200;

% �������: ����չ��
derivatives = [];  % ���ڼ�¼����
ufunc = @fun_0_lorenz;
a = 0;
h = 0.01;
b = length;
n=floor((b-a)/h);       %����
x(1)=a;                 %ʱ�����
y(:,1)=y0;              %����ֵ������������������Ҫע��ά��
for i=1:n               %�����������������ֵ���
    x(i+1)=x(i)+h;
    k1=ufunc(x(i),y(:,i));
    k2=ufunc(x(i)+h/2,y(:,i)+h*k1/2);
    k3=ufunc(x(i)+h/2,y(:,i)+h*k2/2);
    k4=ufunc(x(i)+h,y(:,i)+h*k3);
    y(:,i+1)=y(:,i)+h*(k1+2*k2+2*k3+k4)/6;
    derivatives = [derivatives, (k1+2*k2+2*k3+k4)/6];  % �������Ľ�����м�¼.
end
y = y';


% ��ͼ
figure
plot(y(:,1),y(:,3));  % xz��ͼ�켣
grid on;

% -----------------------------

% ���һ����
the_point = y(end, :);

% ���� ������� ���ú���
[t_n, y_n] = fun_8_RungeKutta(@fun_18_negative_lorenz, 0, 0.01, length, the_point);
y_n = y_n';

% ���� ���ڵ�������
y_n = [the_point];
y_n = y_n';  % תΪʱ������
for i=1:size(derivatives, 2)
    y_n(:, i+1) = y_n(:,i)-h*derivatives(:,end-(i-1));  % (1)Ҫ�����������޸� (2)�����һ���ĵ�����ʼִ��
end
y_n = y_n';  % ת��Ϊʱ������

figure
back_n = size(derivatives, 2);
plot(y_n(1:back_n ,1), y_n(1:back_n ,3));  % xz��ͼ�켣
grid on;

% start point
start_point = y_n(end, :);