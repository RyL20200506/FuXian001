% 官网的线性插值demo, 在复现Yan中, 需要求dot(x), 直接线性插值


x = 0:pi/4:2*pi; 
v = sin(x);

xq = 0:pi/16:2*pi;

figure
vq1 = interp1(x,v,xq);
plot(x,v,'o',xq,vq1,':.');
xlim([0 2*pi]);
title('(Default) Linear Interpolation');


