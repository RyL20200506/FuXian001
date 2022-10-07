% 主要功能与main12, 或者main1差不多, 我们可以先做一个对main1的模仿(因为这还算是最精确的情况)
% 看看能不能用比较少的步数计算

clear
clc

% 获取真实轨迹 - 自定义龙格库塔
Y = [-8, 7, 27, ones(1,26)]; 
end_point = 5000;

% 龙格库塔
n=floor((b-a)/h);       %步数
x(1)=a;                 %时间起点
y(:,1)=y0;              %赋初值，可以是向量，但是要注意维数
for i=1:n               %龙格库塔方法进行数值求解
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

% 画出参数的运动
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




