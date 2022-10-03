% ÂåÂ××È»ìãçº¯Êı
function ret = fun_0_lorenz(t,y)

a=10;
b=8/3;
r=28;

ret = [
    a*(y(2)-y(1));
    r*y(1)-y(2)-y(1)*y(3);
    y(1)*y(2)-b*y(3)];
end