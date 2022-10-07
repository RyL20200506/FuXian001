% ��Ҫ������main12, ����main1���, ���ǿ�������һ����main1��ģ��(��Ϊ�⻹�����ȷ�����)
% �����ܲ����ñȽ��ٵĲ�������

clear
clc

% ��ȡ��ʵ�켣 - �Զ����������
Y = [-8, 7, 27, ones(1,26)]; 
end_point = 5000;

% �������
n=floor((b-a)/h);       %����
x(1)=a;                 %ʱ�����
y(:,1)=y0;              %����ֵ������������������Ҫע��ά��
for i=1:n               %�����������������ֵ���
    x(i+1)=x(i)+h;
    k1=ufunc(x(i),y(:,i), para);
    k2=ufunc(x(i)+h/2,y(:,i)+h*k1/2, para);
    k3=ufunc(x(i)+h/2,y(:,i)+h*k2/2, para);
    k4=ufunc(x(i)+h,y(:,i)+h*k3, para);
    y(:,i+1)=y(:,i)+h*(k1+2*k2+2*k3+k4)/6;
end



[t,y] = fun_8_RungeKutta(@fun_1_lorenz_solver2,0,0.005,end_point,Y);
% [t,y] = fun_8_RungeKutta(@fun_1_lorenz_solver,0,0.005,end_point,Y);
X_n = y';

% �����������˶�
figure
plot(X_n(:,25)) % a
hold on 
plot(X_n(:,26)) % b
hold on 
plot(X_n(:,27)) % r
hold on 
legend('a','b','r') 
% save('mat_1_y', 'y')
para = X_n(end, [25,26,27])

stop




