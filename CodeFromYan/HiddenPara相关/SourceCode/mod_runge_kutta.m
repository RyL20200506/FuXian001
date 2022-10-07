function [t,y]=mod_runge_kutta(ufunc,y0,h,aa,bb,x_ce,x_dot,a,b,r,alfa,learn_rate)
n=floor((bb-aa)/h);       %步数
t(1)=aa;                 %时间起点
y(:,1)=y0;              %赋初值，可以是向量，但是要注意维�?
for i=1:n               %龙格库塔方法进行数�?求解    
    t(i+1)=t(i)+h;    
    k1=ufunc(t(i),y(:,i),x_ce(i,1),x_dot(i,1),a,b,r,alfa,learn_rate);  
    k2=ufunc(t(i)+h/2,y(:,i)+h*k1/2,(x_ce(i,1)+x_ce(i+1,1))/2,(x_dot(i,1)+x_dot(i+1,1))/2,a,b,r,alfa,learn_rate);    
    k3=ufunc(t(i)+h/2,y(:,i)+h*k2/2,(x_ce(i,1)+x_ce(i+1,1))/2,(x_dot(i,1)+x_dot(i+1,1))/2,a,b,r,alfa,learn_rate);   
    k4=ufunc(t(i)+h,y(:,i)+h*k3,x_ce(i+1,1),x_dot(i+1,1),a,b,r,alfa,learn_rate);   
    y(:,i+1)=y(:,i)+h*(k1+2*k2+2*k3+k4)/6;
end