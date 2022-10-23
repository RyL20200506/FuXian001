% ����ѧϰ��Ԫ�������

% ����ѵ����
rng default % for reproducibility
tdata = 0:0.1:10;
ydata = 40*exp(-0.5*tdata) + randn(size(tdata));

% ������Ҫ��С���ĺ���, x��ϣ����ϵĲ���, ����tdata��ȥ, ϣ���õ�ydata
fun = @(x)fun_2_sseval(x,tdata,ydata);

% �������Ų���
x0 = rand(2,1);
bestx = fminsearch(fun,x0);

%����������
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
