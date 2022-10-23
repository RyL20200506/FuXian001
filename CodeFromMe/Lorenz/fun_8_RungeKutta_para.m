% 希望用来实现龙格库塔
% 输入一个函数, 和当前的值, 以及步长, 以及迭代的次数; 


%参数表顺序依次是微分方程组的函数名称，初始值向量，步长，时间起点，时间终点（参数形式参考了ode45函数）
function [time_list,y_list]=fun_8_RungeKutta_para(ufunc,a,h,b,y0, para)
n=floor((b-a)/h);       %步数
time_list(1)=a;                 %时间起点
y_list(:,1)=y0;              %赋初值，可以是向量，但是要注意维数

for i=1:n               %龙格库塔方法进行数值求解
    time_list(i+1)=time_list(i)+h;
    k1=ufunc(time_list(i), y_list(:,i), para);
    k2=ufunc(time_list(i)+h/2, y_list(:,i)+h*k1/2, para);
    k3=ufunc(time_list(i)+h/2, y_list(:,i)+h*k2/2, para);
    k4=ufunc(time_list(i)+h, y_list(:,i)+h*k3, para);
    y_list(:,i+1)=y_list(:,i)+h*(k1+2*k2+2*k3+k4)/6;
end