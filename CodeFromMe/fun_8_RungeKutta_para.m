% 希望用来实现龙格库塔
% 输入一个函数, 和当前的值, 以及步长, 以及迭代的次数; 


%参数表顺序依次是微分方程组的函数名称，初始值向量，步长，时间起点，时间终点（参数形式参考了ode45函数）
function [x,y]=fun_8_RungeKutta_para(ufunc,a,h,b,y0, para)
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