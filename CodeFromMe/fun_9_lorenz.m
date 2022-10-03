% ÂåÂ××È»ìãçº¯Êý
function ret = fun_9_lorenz(t,y, F)

a=10;
b=8/3;
r=28;

a*(y(2)-y(1));
% dot_x = F(t)

ret = [ 
    a*(y(2)-y(1));
    r*y(1)-y(2)-y(1)*y(3);
    y(1)*y(2)-b*y(3)];
end