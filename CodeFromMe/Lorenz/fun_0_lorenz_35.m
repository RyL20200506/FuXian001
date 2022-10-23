% ÂåÂ××È»ìãçº¯Êı
function ret = fun_0_lorenz_35(t,y)

a=15;
b=3;
r=35;

ret = [ 
    a*(y(2)-y(1));
    r*y(1)-y(2)-y(1)*y(3);
    y(1)*y(2)-b*y(3)];
end