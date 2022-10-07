% 结果: 
% (1)代码是没有问题的, 其预言的r值不会出现偏差
% (2)所写的方程跟我写的不一样
% todo: 对比他们之间的区别!

function dy=mod_lor_fun(t,y,x,x1_dot,a,b,r,alfa,learn_rate)
global VAL;
theta=1;
beta=alfa;
beta1=alfa;
betaA=learn_rate;
beta3=learn_rate;
beta4=learn_rate;
beta30=learn_rate;
% beta2=20;
x_ce=x;  % 真实的x
x_dot=x1_dot;  % 真实的x_dot

% f_deta=-2*a*(x_dot-a*(y(1)-x_ce))*theta;

% VAL=[VAL;f_deta];
% z_deta=2*x_ce*f_deta*0.01;
% z_deta=0;
% a1=0.01*((y(1)-x_ce)*(x_dot-y(3)*(y(1)-x_ce)));
a=y(14);
r=y(16);
b=y(18);
y1=y(20);
z1=y(22);
% learn=log(abs(x_dot-a*(y(11)-x_ce)));
% if learn>0
% beta3=0.01*learn;
% beta4=0.01*learn;
% else
%    beta3=0.01;
%    beta4=0.01; 
% end

% a=y(13);
% r=y(14);
% b=y(15);
% y1=y(16);
% z1=y(17);

gam=1;

dy=[
    % 可能是对a b r e1 e2的求导?
    -y(1)-x_ce*y(2)-beta*y(1);  % y(1): 预测的y对a的导数
    x_ce*y(1)-b*y(2)-beta*y(2);  % y(2): 预测的z对a的导数
    
    x_ce-y(3)-x_ce*y(4)-beta*y(3);  % y(3): 预测的y对r的导数
    x_ce*y(3)-b*y(4)-beta*y(4);  % y(4):  预测的z对r的导数
    
    -y(5)-x_ce*y(6)-beta*y(5);  % y(5): 预测的y对b的导数
    x_ce*y(5)-y(12)-b*y(6)-beta*y(6);  % y(6):  预测的z对b的导数
    
    -y(7)-x_ce*y(8)+1-beta*y(7);  % y(7): 预测的y对e1的导数!
    x_ce*y(7)-b*y(8)-beta*y(8);  % y(8); 预测的z对e1的导数?
    
    -y(9)-x_ce*y(10)-beta*y(9);  % y(9): 预测的y对e2的导数?
    x_ce*y(9)-b*y(10)+1-beta*y(10);  % y(10): 预测的z对e2的导数!
    
    % 预测的y和z的动力学
    r*x_ce-y(11)-x_ce*y(12)+y1;  % y(11): 表示预测的y: y1对应的是e1, 因为这就是拟合系统的的动力学方程.
    x_ce*y(11)-b*y(12)+z1;  % y(12): 表示预测的z, z1对应的是e2
    
    % 应该是Delta对其他参数的求导
    -1*beta*y(13) + 2*(y(11)-x_ce) * (x_dot-a*(y(11)-x_ce)-a^2*y(1)) + 2*a*x_dot*y(1);  % y(13): hat_a关于Delta的梯度
    betaA*y(13);  % y(14);  应该就是梯度下降法: 表示预测的a的动力学
    
    -1*beta*y(15) + 2*a*(x_dot-a*(y(11)-x_ce))*y(3);  % y(15)
    beta3*y(15);  % y(16);  表示预测的r的动力学
    
    -1*beta*y(17) + 2*a*(x_dot-a*(y(11)-x_ce))*y(5);  % y(17)
    beta4*y(17);  % y(18);  表示预测的b的动力学
    
    -1*beta1*y(19)+2*a*(x_dot-a*(y(11)-x_ce))*y(7)-2*gam*y1;  % y(19): delta 对y1求导, 即对e1求导
    beta30*y(19);  % y(20); y1的动力学, 即e1的动力学 
    
    -1*beta*y(21)+2*a*(x_dot-a*(y(11)-x_ce))*y(9)-2*gam*z1;  % y(21): 需要y对e2求导 其中用到了y(9) 那应该是
    beta30*y(21);    % y(22):  z1的动力学, 即e2的动力学
    
];


% dy=[
%     -y(1)-x_ce*y(1)-beta*y(1);
%     x_ce*y(1)-b*y(2)-beta*y(2);
%     
%     x_ce-y(3)-x_ce*y(4)-beta*y(3);
%     x_ce*y(3)-b*y(4)-beta*y(4);
%     
%     -y(5)-x_ce*y(6)-beta*y(5);
%     x_ce*y(5)-y(12)-b*y(6)-beta*y(6);
%     
%     -y(7)-x_ce*y(8)+1-beta*y(7);
%     x_ce*y(7)-b*y(8)-beta*y(8);
%     
%     -y(9)-x_ce*y(10)-beta*y(9);
%     x_ce*y(9)-b*y(10)+1-beta*y(10);
%     
%     r*x_ce-y(11)-x_ce*y(12)+y1;
%     x_ce*y(11)-b*y(12)+z1;
%     
%    -1*beta*y(13)+2*(y(11)-x_ce)*(x_dot-a*(y(11)-x_ce)-a^2*y(1))+2*a*x_dot*y(1);
%     
%     
%    -1*beta*y(14)+2*a*(x_dot-a*(y(11)-x_ce))*y(3);
%   
%     
%    -1*beta*y(14)+2*a*(x_dot-a*(y(11)-x_ce))*y(5);
%    
%     
%    -1*beta*y(14)+2*a*(x_dot-a*(y(11)-x_ce))*y(7)-2*gam*y1;
%   
%     
%    -1*beta*y(14)+2*a*(x_dot-a*(y(11)-x_ce))*y(9)-2*gam*z1;
%  
%     
% ];





end