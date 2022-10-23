% 洛伦兹混沌 求解参数的动力学
% fun_9: 与之前的区别在于, 
function ret = fun_lorenz_solver_direction_positive(t,Y, FXdot, FX)

% 参数
direction = 1;  % 表示方向朝着正向流逝
gamma=0.0035;
alpha_m=3;
beta=1;

% 命名
haty = Y(1);
hatz = Y(2);

haty_hata = Y(3);
haty_hatb = Y(4);
haty_hatr = Y(5);
haty_e1 = Y(6);
haty_e2 = Y(7);

hatz_hata = Y(8);
hatz_hatb = Y(9);
hatz_hatr = Y(10);
hatz_e1 = Y(11);
hatz_e2 = Y(12);

D_hata = Y(13);
D_hatb = Y(14);
D_hatr = Y(15);
D_e1 = Y(16);
D_e2 = Y(17);

hata = Y(18);
hatb = Y(19);
hatr = Y(20);
e1 = Y(21);
e2 = Y(22);

% 真实轨迹数据
x = FX(t);  % 这里我们带入真实的X
dot_x = FXdot(t);  % 或者, 用插值的函数作为为最后的结果;(显然如果我们用多项式2次插值, 应该可以得到更好的效果) 

ret=[
% % % 拟合系统 无误
direction * ( hatr*x - haty - x*hatz + e1 );  % haty
direction * ( x*haty - hatb*hatz + e2 );  % hatz

% haty的偏导
direction * ( -haty_hata + (-x)*hatz_hata );  % haty_hata
direction * ( -haty_hatb +(-x)*hatz_hatb ); 
direction * ( x - haty_hatr + (-x)*hatz_hatr );  
direction * ( -haty_e1 + (-x)*hatz_e1 + 1 );  
direction * ( -haty_e2 + (-x)*hatz_e2 );  

% hatz的偏导
direction * ( x*haty_hata - hatb*hatz_hata );  % Y(15):Y(15)
direction * ( x*haty_hatb + (-1)*(hatz + hatb*hatz_hatb) );  % Y(16):Y(16)
direction * ( x*haty_hatr - hatb*hatz_hatr );  % Y(17):Y(17)
direction * ( x*haty_e1 - hatb*hatz_e1 );  % Y(6)_hat{Y(28)}:Y(6)_hat{Y(28)}
direction * ( x*haty_e2 - hatb*hatz_e2 + 1 );  % Y(6)_hat{Y(29)}:Y(6)_hat{Y(29)}

% Delta的导数
-alpha_m*D_hata + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ( (haty - x) + hata * haty_hata ) ;
-alpha_m*D_hatb + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ( hata * haty_hatb );
-alpha_m*D_hatr + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ( hata * haty_hatr );
-alpha_m*D_e1 + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ( hata * haty_e1 ) + 2*beta* e1 ;  
-alpha_m*D_e2 + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ( hata * haty_e2 ) + 2*beta* e2 ;  

% 参数动力学 应该无误
-2*gamma * D_hata;  % Y(25): Y(25)
-2*gamma * D_hatb;  % Y(26): Y(26)
-2*gamma * D_hatr;  % Y(27): Y(27)
-2*gamma * D_e1;  % Y(28): hat{Y(28)}
-2*gamma * D_e2;  % Y(29): hat{Y(29)}


%-----------
% direction * ( hatr*x - haty - x*hatz + e1 );  % haty
% direction * ( x*haty - hatb*hatz + e2 );  % hatz
% 
% % haty的偏导
% direction * ( -haty_hata + (-x)*hatz_hata );  % haty_hata
% direction * ( -haty_hatb +(-x)*hatz_hatb ); 
% direction * ( x - haty_hatr + (-x)*hatz_hatr );  
% direction * ( -haty_e1 + (-x)*hatz_e1 + 1 );  
% direction * ( -haty_e2 + (-x)*hatz_e2 );  
% 
% % hatz的偏导
% direction * ( x*haty_hata - hatb*hatz_hata );  % Y(15):Y(15)
% direction * ( x*haty_hatb + (-1)*(hatz + hatb*hatz_hatb) );  % Y(16):Y(16)
% direction * ( x*haty_hatr - hatb*hatz_hatr );  % Y(17):Y(17)
% direction * ( x*haty_e1 - hatb*hatz_e1 );  % Y(6)_hat{Y(28)}:Y(6)_hat{Y(28)}
% direction * ( x*haty_e2 - hatb*hatz_e2 + 1 );  % Y(6)_hat{Y(29)}:Y(6)_hat{Y(29)}
% 
% % Delta的导数
% -alpha_m*D_hata + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ((haty - x) + 2*hata * haty_hata ) ;
% -alpha_m*D_hatb + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ( 2*hata * haty_hatb );
% -alpha_m*D_hatr + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ( 2*hata * haty_hatr );
% -alpha_m*D_e1 + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction *  (2*hata * haty_e1 ) + 2*beta* e1 ;  
% -alpha_m*D_e2 + (-2)*( dot_x - direction*( hata*(haty-x)) ) * direction * ( 2*hata * haty_e2 ) + 2*beta* e2 ;  
% 
% % 参数动力学 应该无误
% -2*gamma * D_hata;  % Y(25): Y(25)
% -2*gamma * D_hatb;  % Y(26): Y(26)
% -2*gamma * D_hatr;  % Y(27): Y(27)
% -2*gamma * D_e1;  % Y(28): hat{Y(28)}
% -2*gamma * D_e2;  % Y(29): hat{Y(29)}

];

end