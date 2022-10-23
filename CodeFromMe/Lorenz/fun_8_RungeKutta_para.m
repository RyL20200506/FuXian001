% ϣ������ʵ���������
% ����һ������, �͵�ǰ��ֵ, �Լ�����, �Լ������Ĵ���; 


%������˳��������΢�ַ�����ĺ������ƣ���ʼֵ������������ʱ����㣬ʱ���յ㣨������ʽ�ο���ode45������
function [time_list,y_list]=fun_8_RungeKutta_para(ufunc,a,h,b,y0, para)
n=floor((b-a)/h);       %����
time_list(1)=a;                 %ʱ�����
y_list(:,1)=y0;              %����ֵ������������������Ҫע��ά��

for i=1:n               %�����������������ֵ���
    time_list(i+1)=time_list(i)+h;
    k1=ufunc(time_list(i), y_list(:,i), para);
    k2=ufunc(time_list(i)+h/2, y_list(:,i)+h*k1/2, para);
    k3=ufunc(time_list(i)+h/2, y_list(:,i)+h*k2/2, para);
    k4=ufunc(time_list(i)+h, y_list(:,i)+h*k3, para);
    y_list(:,i+1)=y_list(:,i)+h*(k1+2*k2+2*k3+k4)/6;
end