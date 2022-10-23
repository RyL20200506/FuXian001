% ÂåÂ××È»ìãçº¯Êı: 
function ret = fun_0_lorenz_para(t,y, para)

a=para(1);
b=para(2);
r=para(3);

ret = [ 
    a*(y(2)-y(1));
    r*y(1)-y(2)-y(1)*y(3);
    y(1)*y(2)-b*y(3)];
end