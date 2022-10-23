% 洛伦兹混沌函数: 带有噪声的版本
function ret = fun_11_lorenz_noise(t,y)

a=10;
b=8/3;
r=28;
D=0;  % 噪声强度
ret = [ 
    a*(y(2)-y(1)) + D*randn(1);
    r*y(1)-y(2)-y(1)*y(3) + D*randn(1);
    y(1)*y(2)-b*y(3)] + D*randn(1);
end