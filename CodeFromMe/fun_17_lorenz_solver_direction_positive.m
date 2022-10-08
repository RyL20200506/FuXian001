% 洛伦兹混沌 求解参数的动力学
% fun_9: 与之前的区别在于, 
function ret = fun_17_lorenz_solver_direction_positive(t,Y, FXdot, FX)

direction = 1;  % 表示方向朝着正向流逝

% 根据t定位dot_X
if (mod(t,1)==0)
    t  % 用于验证时间调用该函数时所是用的步长, 发现并不是等步长的
end

% 从Y里面获取值
% 赋值
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

% 下面是一些参数
% a=10;
% b=8/3;
% r=28;
gamma=0.0015;
alpha=3;
beta=2;

% dot_x = a*(Y(2)-Y(1));  % 直接用真实的结果作为dot_x
dot_x = FXdot(t);  % 或者, 用插值的函数作为为最后的结果;(显然如果我们用多项式2次插值, 应该可以得到更好的效果) 
x = FX(t);  % 这里我们带入真实的X

ret=[

% 原系统的方程: 其实不需要原方程;应该也依然work
0; % a*(Y(2)-Y(1));  % Y(1):Y(1), a*(Y(2)-Y(1))
0; % r*Y(1)-Y(2)-Y(1)*Y(3);  % Y(2):Y(2), r*Y(1)-Y(2)-Y(1)*Y(3)
0; % Y(1)*Y(2)-b*Y(3);  % Y(3):Y(3), Y(1)*Y(2)-b*Y(3)

% 拟合系统
direction * ( hata*(haty-x) ) ;  % Y(4):Y(4), Y(25)*(Y(5)-Y(1))
direction * ( hatr*x - haty - x*hatz + e1 );  % Y(5):Y(5), 
direction * ( x*haty - hatb*hatz + e2 );  % Y(6):Y(6)

% hatx的偏导
direction * ( (haty - x) + hata * haty_hata ); % Y(7):Y(7)
direction * ( hata * haty_hatb );  % Y(8):Y(4)_b
direction * ( hata * haty_hatr );  % Y(9):Y(4)_r = Y(25) * hat(Y(2))_r
direction * ( hata * haty_e1 );
direction * ( hata * haty_e2 );

% haty的偏导
direction * ( -haty_hata + (-x)*hatz_hata );  % Y(10):Y(10) = -Y(5)_a
direction * ( -haty_hatb +(-x)*hatz_hatb );  % Y(11):Y(11)
direction * ( x - haty_hatr + (-x)*hatz_hatr );  % Y(12):Y(12)
direction * ( -haty_e1 + (-x)*hatz_e1 + 1 );  % Y(5)_hat{Y(28)}:Y(5)_hat{Y(28)}
direction * ( -haty_e2 + (-x)*hatz_e2 );  % Y(5)_hat{Y(29)}:Y(5)_hat{Y(29)}

% hatz的偏导
direction * ( x*haty_hata - hatb*hatz_hata );  % Y(15):Y(15)
direction * ( x*haty_hatb + (-1)*(hatz + hatb*hatz_hatb) );  % Y(16):Y(16)
direction * ( x*haty_hatr - hatb*hatz_hatr );  % Y(17):Y(17)
direction * ( x*haty_e1 - hatb*hatz_e1 );  % Y(6)_hat{Y(28)}:Y(6)_hat{Y(28)}
direction * ( x*haty_e2 - hatb*hatz_e2 + 1 );  % Y(6)_hat{Y(29)}:Y(6)_hat{Y(29)}

% Delta的导数
-alpha*D_hata + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ( (haty - x) + hata * haty_hata ) ;
-alpha*D_hatb + (-2)*( dot_x - direction*( hata*(haty-x)) )* direction * ( hata * haty_hatb );
-alpha*D_hatr + (-2)*( dot_x - direction*( hata*(haty-x)) )* direction * ( hata * haty_hatr );
-alpha*D_e1 + (-2)*( dot_x - direction*( hata*(haty-x)) )* direction * ( hata * haty_e1 ) + 2*beta* e1 ;  
-alpha*D_e2 + (-2)*( dot_x - direction*( hata*(haty-x)) )* direction * ( hata * haty_e2 ) + 2*beta* e2 ;  

% 参数动力学 应该无误
-2*gamma * D_hata;  % Y(25): Y(25)
-2*gamma * D_hatb;  % Y(26): Y(26)
-2*gamma * D_hatr;  % Y(27): Y(27)
-2*gamma * D_e1;  % Y(28): hat{Y(28)}
-2*gamma * D_e2;  % Y(29): hat{Y(29)}
];

end