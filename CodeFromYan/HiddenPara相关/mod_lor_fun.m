% ���: 
% (1)������û�������, ��Ԥ�Ե�rֵ�������ƫ��
% (2)��д�ķ��̸���д�Ĳ�һ��
% todo: �Ա�����֮�������!

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
x_ce=x;  % ��ʵ��x
x_dot=x1_dot;  % ��ʵ��x_dot

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
    % �����Ƕ�a b r e1 e2����?
    -y(1)-x_ce*y(2)-beta*y(1);  % y(1): Ԥ���y��a�ĵ���
    x_ce*y(1)-b*y(2)-beta*y(2);  % y(2): Ԥ���z��a�ĵ���
    
    x_ce-y(3)-x_ce*y(4)-beta*y(3);  % y(3): Ԥ���y��r�ĵ���
    x_ce*y(3)-b*y(4)-beta*y(4);  % y(4):  Ԥ���z��r�ĵ���
    
    -y(5)-x_ce*y(6)-beta*y(5);  % y(5): Ԥ���y��b�ĵ���
    x_ce*y(5)-y(12)-b*y(6)-beta*y(6);  % y(6):  Ԥ���z��b�ĵ���
    
    -y(7)-x_ce*y(8)+1-beta*y(7);  % y(7): Ԥ���y��e1�ĵ���!
    x_ce*y(7)-b*y(8)-beta*y(8);  % y(8); Ԥ���z��e1�ĵ���?
    
    -y(9)-x_ce*y(10)-beta*y(9);  % y(9): Ԥ���y��e2�ĵ���?
    x_ce*y(9)-b*y(10)+1-beta*y(10);  % y(10): Ԥ���z��e2�ĵ���!
    
    % Ԥ���y��z�Ķ���ѧ
    r*x_ce-y(11)-x_ce*y(12)+y1;  % y(11): ��ʾԤ���y: y1��Ӧ����e1, ��Ϊ��������ϵͳ�ĵĶ���ѧ����.
    x_ce*y(11)-b*y(12)+z1;  % y(12): ��ʾԤ���z, z1��Ӧ����e2
    
    % Ӧ����Delta��������������
    -1*beta*y(13) + 2*(y(11)-x_ce) * (x_dot-a*(y(11)-x_ce)-a^2*y(1)) + 2*a*x_dot*y(1);  % y(13): hat_a����Delta���ݶ�
    betaA*y(13);  % y(14);  Ӧ�þ����ݶ��½���: ��ʾԤ���a�Ķ���ѧ
    
    -1*beta*y(15) + 2*a*(x_dot-a*(y(11)-x_ce))*y(3);  % y(15)
    beta3*y(15);  % y(16);  ��ʾԤ���r�Ķ���ѧ
    
    -1*beta*y(17) + 2*a*(x_dot-a*(y(11)-x_ce))*y(5);  % y(17)
    beta4*y(17);  % y(18);  ��ʾԤ���b�Ķ���ѧ
    
    -1*beta1*y(19)+2*a*(x_dot-a*(y(11)-x_ce))*y(7)-2*gam*y1;  % y(19): delta ��y1��, ����e1��
    beta30*y(19);  % y(20); y1�Ķ���ѧ, ��e1�Ķ���ѧ 
    
    -1*beta*y(21)+2*a*(x_dot-a*(y(11)-x_ce))*y(9)-2*gam*z1;  % y(21): ��Ҫy��e2�� �����õ���y(9) ��Ӧ����
    beta30*y(21);    % y(22):  z1�Ķ���ѧ, ��e2�Ķ���ѧ
    
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