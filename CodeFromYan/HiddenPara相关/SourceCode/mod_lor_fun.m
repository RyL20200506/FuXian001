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
x_ce=x;
x_dot=x1_dot;

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
    -y(1)-x_ce*y(2)-beta*y(1);
    x_ce*y(1)-b*y(2)-beta*y(2);
    
    x_ce-y(3)-x_ce*y(4)-beta*y(3);
    x_ce*y(3)-b*y(4)-beta*y(4);
    
    -y(5)-x_ce*y(6)-beta*y(5);
    x_ce*y(5)-y(12)-b*y(6)-beta*y(6);
    
    -y(7)-x_ce*y(8)+1-beta*y(7);
    x_ce*y(7)-b*y(8)-beta*y(8);
    
    -y(9)-x_ce*y(10)-beta*y(9);
    x_ce*y(9)-b*y(10)+1-beta*y(10);
    
    r*x_ce-y(11)-x_ce*y(12)+y1;
    x_ce*y(11)-b*y(12)+z1;
    
    -1*beta*y(13)+2*(y(11)-x_ce)*(x_dot-a*(y(11)-x_ce)-a^2*y(1))+2*a*x_dot*y(1);
    betaA*y(13);
    
    -1*beta*y(15)+2*a*(x_dot-a*(y(11)-x_ce))*y(3);
    beta3*y(15);
    
    -1*beta*y(17)+2*a*(x_dot-a*(y(11)-x_ce))*y(5);
    beta4*y(17);
    
    -1*beta1*y(19)+2*a*(x_dot-a*(y(11)-x_ce))*y(7)-2*gam*y1;
    beta30*y(19);
    
    -1*beta*y(21)+2*a*(x_dot-a*(y(11)-x_ce))*y(9)-2*gam*z1;
    beta30*y(21);
    
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