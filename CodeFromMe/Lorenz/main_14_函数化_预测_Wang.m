%% HOPE
% ��(1)��֪��������ʽ, (2)ֻ����һ�����ݵ������, (3)�����Ƚ�׼ȷ��Ԥ��
close all;clear;clc;

% ���ɴ�Ԥ������
[t1,h1]=fun_8_RungeKutta(@fun_0_lorenz, 0, 0.01, 200, [12 4 0]);
plot3(h1(1,:),h1(2,:),h1(3,:),'r'); grid on;
history = h1';
save('mat_14_history.mat', 'history')

% ׼������
history_train = history(1:10000, :);
predict_point = history_train(end, :);
predict_length = 5000;

% ��ʼԤ��
future_values = fun_14_GeneralPredict_Wang(history_train, predict_point, predict_length);

% ��Ԥ���ͼ: Ч��һ��, ���Ԥ��
history_test = history(10000:end, :)
plot(history_test(:,1))
hold on
plot(future_values(:,1))

stop;

%% ��������: ��¼A3��������
close all;clear;clc;
delta=10;b=8/3;r=28;
N=15000;
% x0=[-1;3;4];
x0=[-8;7;27];
X_n=[x0];
t0=0.01;%��λ����ʱ��
m=1;
D=0;%��������
% D=0;
b1=sqrt(2*D*t0);b2=1/2*b1*t0;b3=1/sqrt(3)*b2;
w1=randn(N,3);w2=randn(N,3);
for i=1:N
    s1=b1*w1(i,:);
    s2=b2*w1(i,:)+b3*w2(i,:);
    df2=[0;(r-X_n(3,i)-1-X_n(1,i));(X_n(2,i)+X_n(1,i)-b)];
    fx=delta*(X_n(2,i)-X_n(1,i));
    fy=r*X_n(1,i)-X_n(2,i)-X_n(1,i)*X_n(3,i);
    fz=X_n(1,i)*X_n(2,i)-b*X_n(3,i);
    df1=[-delta*fx+delta*fy;(r-X_n(3,i))*fx-fy-X_n(1,i)*fz;X_n(2,i)*fx+X_n(1,i)*fy-b*fz]; 
    X_t=X_n(:,i)+[fx;fy;fz]*t0+1/2*t0*t0*df1+df2.*s2'+s1';
    X_n=[X_n,X_t];
end
history = X_n';

% ׼������
history_train = history(1:10000, :);
predict_point = history_train(end, :);
predict_length = 4000;

% ��ʼԤ��
future_values = fun_14_GeneralPredict_Wang(history_train, predict_point, predict_length);

% ��Ԥ���ͼ: Ԥ�ⷶΧ2800
figure
history_test = history(10000:end, :)
plot(history_test(:,1))
hold on
plot(future_values(:,1))

stop;

%% ����Ϊxie�����ݼ�
close all;clear;clc;

% ���ɴ�Ԥ������
load mat_14_history_from_xie.mat;

% ׼������
history_train = history(1:8000, :);
start_from = 1;
predict_point = history_train(start_from, :);
predict_length = 2000;

% ��ʼԤ��
future_values = fun_14_GeneralPredict_Wang(history_train, predict_point, predict_length);

% ��Ԥ���ͼ
figure
history_test = history(start_from:end, :)
plot(history_test(:,1))
hold on
plot(future_values(:,1))

stop;








