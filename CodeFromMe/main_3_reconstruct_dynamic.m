% ����: Wang - �ع�����ѧ
% Wang ����

close all;clear;clc;

% ����ѵ����
y0 = [1,1,1];  % ��ʼֵ
[t, real_trajectory] = ode45('fun_0_lorenz',[0,100],y0); % ����֮ǰ�ĺ��������ʵ�켣
X = real_trajectory(1:end-1, :);  % �ָ������
Y = real_trajectory(2:end, :);  % �ٷָ����Ԥ���ǩ

% ������Ҫ��С���ĺ���, x��ϣ����ϵĲ���, ����tdata��ȥ, ϣ���õ�ydata
fun = @(PARA)fun_3_costfun(PARA,X,Y);

% �������Ų���
PARA_0 = rand(48,1);  % ��ʼ����
options = optimset('MaxFunEvals',10000000);
best_PARA = fminsearch(fun,PARA_0, options);
% ���µ����˼���
best_PARA = fminsearch(fun,best_PARA, options);  
best_PARA = fminsearch(fun,best_PARA, options); 
best_PARA = fminsearch(fun,best_PARA, options); 
best_PARA = fminsearch(fun,best_PARA, options); 

% --------------------------------------------------
% ����ѵ����
% ����ѵ����
y0 = [1,1,1];  % ��ʼֵ
[t, real_trajectory] = ode45('fun_0_lorenz',[20,50],y0); % ����֮ǰ�ĺ��������ʵ�켣
X = real_trajectory(1:end-1, :);  % �ָ������
Y = real_trajectory(2:end, :);  % �ٷָ����Ԥ���ǩ

% ������Ҫ��С���ĺ���, x��ϣ����ϵĲ���, ����tdata��ȥ, ϣ���õ�ydata
fun = @(PARA)fun_3_costfun(PARA,X,Y);

% �������Ų���
PARA_0 = rand(48,1);  % ��ʼ����
options = optimset('MaxFunEvals',10000000);
% ���µ�������
best_PARA = fminsearch(fun,best_PARA, options);  
best_PARA = fminsearch(fun,best_PARA, options); 


% --------------------------------------------------
% --------------------------------------------------
% ����������: �Ϳ���ģ�͵ľ��� 
X0 = X(1,:);
predict_result = fun_3_4_predict(best_PARA, X0, 1000); 

% 
plot(predict_result(:, 1), '*');
hold on
plot(X(1:1000, 1),'r');

% ����best para�����Ч��
A = bestx(1);
lambda = bestx(2);
yfit = A*exp(-lambda*tdata);
plot(tdata,ydata,'*');
hold on
plot(tdata,yfit,'r');
xlabel('tdata')
ylabel('Response Data and Curve')
title('Data and Best Fitting Exponential Curve')
legend('Data','Fitted Curve')
hold off
