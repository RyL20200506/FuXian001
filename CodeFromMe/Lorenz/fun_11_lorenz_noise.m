% �����Ȼ��纯��: ���������İ汾
function ret = fun_11_lorenz_noise(t,y)

a=10;
b=8/3;
r=28;
D=0;  % ����ǿ��
ret = [ 
    a*(y(2)-y(1)) + D*randn(1);
    r*y(1)-y(2)-y(1)*y(3) + D*randn(1);
    y(1)*y(2)-b*y(3)] + D*randn(1);
end