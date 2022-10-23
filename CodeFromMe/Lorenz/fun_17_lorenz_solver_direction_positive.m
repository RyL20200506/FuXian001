% �����Ȼ��� �������Ķ���ѧ
% fun_9: ��֮ǰ����������, 
function ret = fun_17_lorenz_solver_direction_positive(t,Y, FXdot, FX)

direction = 1;  % ��ʾ��������������

% ����t��λdot_X
if (mod(t,1)==0)
    t  % ������֤ʱ����øú���ʱ�����õĲ���, ���ֲ����ǵȲ�����
end

% ��Y�����ȡֵ
% ��ֵ
x = Y(1);
y = Y(2);
z = Y(3);
hatx = Y(4);
haty = Y(5);
hatz = Y(6);

hatx_hata = Y(7);
hatx_hatb = Y(8);
hatx_hatr = Y(9);
hatx_e1 = Y(10);
hatx_e2 = Y(11);

haty_hata = Y(12);
haty_hatb = Y(13);
haty_hatr = Y(14);
haty_e1 = Y(15);
haty_e2 = Y(16);

hatz_hata = Y(17);
hatz_hatb = Y(18);
hatz_hatr = Y(19);
hatz_e1 = Y(20);
hatz_e2 = Y(21);

D_hata = Y(22);
D_hatb = Y(23);
D_hatr = Y(24);
D_e1 = Y(25);
D_e2 = Y(26);

hata = Y(27);
hatb = Y(28);
hatr = Y(29);
e1 = Y(30);
e2 = Y(31);

% ������һЩ����
% a=10;
% b=8/3;
% r=28;
gamma=0.0015;
alpha_m=3;
beta=2;

% dot_x = a*(Y(2)-Y(1));  % ֱ������ʵ�Ľ����Ϊdot_x
dot_x = FXdot(t);  % ����, �ò�ֵ�ĺ�����ΪΪ���Ľ��;(��Ȼ��������ö���ʽ2�β�ֵ, Ӧ�ÿ��Եõ����õ�Ч��) 
x = FX(t);  % �������Ǵ�����ʵ��X

ret=[

% ԭϵͳ�ķ���: ��ʵ����Ҫԭ����;Ӧ��Ҳ��Ȼwork
0; % a*(Y(2)-Y(1));  % Y(1):Y(1), a*(Y(2)-Y(1))
0; % r*Y(1)-Y(2)-Y(1)*Y(3);  % Y(2):Y(2), r*Y(1)-Y(2)-Y(1)*Y(3)
0; % Y(1)*Y(2)-b*Y(3);  % Y(3):Y(3), Y(1)*Y(2)-b*Y(3)

% ���ϵͳ ����
direction * ( hata*(haty-x) ) ;  % hatx
direction * ( hatr*x - haty - x*hatz + e1 );  % haty
direction * ( x*haty - hatb*hatz + e2 );  % hatz

% hatx��ƫ�� ����
direction * ( (haty - x) + hata * haty_hata ); % hatx_hata
direction * ( hata * haty_hatb );  % hatx_hatb
direction * ( hata * haty_hatr );  % hatx_hatr
direction * ( hata * haty_e1 );
direction * ( hata * haty_e2 );

% haty��ƫ��
direction * ( -haty_hata + (-x)*hatz_hata );  % haty_hata
direction * ( -haty_hatb +(-x)*hatz_hatb ); 
direction * ( x - haty_hatr + (-x)*hatz_hatr );  
direction * ( -haty_e1 + (-x)*hatz_e1 + 1 );  
direction * ( -haty_e2 + (-x)*hatz_e2 );  

% hatz��ƫ��
direction * ( x*haty_hata - hatb*hatz_hata );  % Y(15):Y(15)
direction * ( x*haty_hatb + (-1)*(hatz + hatb*hatz_hatb) );  % Y(16):Y(16)
direction * ( x*haty_hatr - hatb*hatz_hatr );  % Y(17):Y(17)
direction * ( x*haty_e1 - hatb*hatz_e1 );  % Y(6)_hat{Y(28)}:Y(6)_hat{Y(28)}
direction * ( x*haty_e2 - hatb*hatz_e2 + 1 );  % Y(6)_hat{Y(29)}:Y(6)_hat{Y(29)}

% Delta�ĵ���
-alpha_m*D_hata + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ((haty - x) + 2*hata * haty_hata ) ;
-alpha_m*D_hatb + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * (2* hata * haty_hatb );
-alpha_m*D_hatr + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * (2* hata * haty_hatr );
-alpha_m*D_e1 + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ( hata * haty_e1 ) + 2*beta* e1 ;  
-alpha_m*D_e2 + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ( hata * haty_e2 ) + 2*beta* e2 ;  

% ��������ѧ Ӧ������
-2*gamma * D_hata;  % Y(25): Y(25)
-2*gamma * D_hatb;  % Y(26): Y(26)
-2*gamma * D_hatr;  % Y(27): Y(27)
-2*gamma * D_e1;  % Y(28): hat{Y(28)}
-2*gamma * D_e2;  % Y(29): hat{Y(29)}
];

end