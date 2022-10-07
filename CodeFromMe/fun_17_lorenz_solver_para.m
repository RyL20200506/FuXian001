% 洛伦兹混沌 求解参数的动力学
% fun_9: 与之前的区别在于, 
function ret = fun_17_lorenz_solver_para(t,Y, FXdot, FX)

% x:y(1), y:y(2), z:y(3), 
% hat{x}: y(4), hat{y}: y(5), hat{z}: y(6)
% hat{a}: y(7), hat{b}: y(8), hat{r}: y(9), 
% e_1:y(10), e_2:y(11)
% G_k_a: y(12), G_k_b: y(13), G_k_r: y(14)
% G_k_e_1: y(15), G_k_e_2: y(16)
% H_k_a: y(17), H_k_b: y(18), H_k_r: y(19)
% H_k_e_1: y(20), H_k_e_2: y(21)
% 公式参考:http://127.0.0.1:1111/center/call_exp_general_function/ftn.url_open_view_via_id('20220912-1853-4398-4377-000000000000')

% 根据t定位dot_X
if (mod(t,1)==0)
    t  % 用于验证时间调用该函数时所是用的步长, 发现并不是等步长的
end

% 下面是一些参数
a=10;
b=8/3;
r=28;
gamma=0.0015;
alpha=3;
beta=2;

dot_x = a*(Y(2)-Y(1));  % 直接用真实的结果作为dot_x
dot_x = FXdot(t);  % 或者, 用插值的函数作为为最后的结果;(显然如果我们用多项式2次插值, 应该可以得到更好的效果) 
x = FX(t);  % 这里我们带入真实的X

ret=[

% 原系统的方程: 其实不需要原方程;应该也依然work
0; % a*(Y(2)-Y(1));  % Y(1):Y(1), a*(Y(2)-Y(1))
0; % r*Y(1)-Y(2)-Y(1)*Y(3);  % Y(2):Y(2), r*Y(1)-Y(2)-Y(1)*Y(3)
0; % Y(1)*Y(2)-b*Y(3);  % Y(3):Y(3), Y(1)*Y(2)-b*Y(3)

% 感觉无误 注意这里带入的是真实的Y(1)
Y(25)*(Y(5)-x);  % Y(4):Y(4), Y(25)*(Y(5)-Y(1))
Y(27)*x-Y(5)-x*Y(6)+Y(28);  % Y(5):Y(5), 
x*Y(5)-Y(26)*Y(6)+Y(29);  % Y(6):Y(6)

% 感觉无误
1*(Y(5)-x)+Y(25)*Y(10); % Y(7):Y(7)
Y(25)*Y(11);  % Y(8):Y(4)_b
Y(25)*Y(12);  % Y(9):Y(4)_r = Y(25) * hat(Y(2))_r

% 感觉无误
-Y(10)+(-x)*Y(15);  % Y(10):Y(10) = -Y(5)_a
-Y(11)+(-x)*Y(16);  % Y(11):Y(11)
x-Y(12)+(-x)*Y(17);  % Y(12):Y(12)
-Y(13)+(-x)*Y(18)+1;  % Y(5)_hat{Y(28)}:Y(5)_hat{Y(28)}
-Y(14)+(-x)*Y(19);  % Y(5)_hat{Y(29)}:Y(5)_hat{Y(29)}

% 感觉无误
x*Y(10)-Y(26)*Y(15);  % Y(15):Y(15)
x*Y(11)+(-1)*(Y(6)+Y(26)*Y(16));  % Y(16):Y(16)
x*Y(12)-Y(26)*Y(17);  % Y(17):Y(17)
x*Y(13)-Y(26)*Y(18);  % Y(6)_hat{Y(28)}:Y(6)_hat{Y(28)}
x*Y(14)-Y(26)*Y(19)+1;  % Y(6)_hat{Y(29)}:Y(6)_hat{Y(29)}

% 参数的动力学 感觉无误
-alpha*Y(20)+(-2)*( dot_x - (Y(25)*(Y(5)-x)) )*(2*Y(25)*Y(10)+(Y(5)-x));
-alpha*Y(21)+(-2)*( dot_x - (Y(25)*(Y(5)-x)) )*(2*Y(25)*Y(11));
-alpha*Y(22)+(-2)*( dot_x - (Y(25)*(Y(5)-x)) )*(2*Y(25)*Y(12));

% 误差梯度动力学 感觉无误
-alpha*Y(23)+(-2)*( dot_x - (Y(25)*(Y(5)-x)) )*Y(25)*Y(13)+2*beta* Y(28); % 这样就比较准了! 但是我还没想明白为什么要这样更改....
-alpha*Y(24)+(-2)*( dot_x - (Y(25)*(Y(5)-x)) )*Y(25)*Y(14)+2*beta* Y(29);

% 参数动力学 应该无误
-2*gamma * Y(20);  % Y(25): Y(25)
-2*gamma * Y(21);  % Y(26): Y(26)
-2*gamma * Y(22);  % Y(27): Y(27)
-2*gamma * Y(23);  % Y(28): hat{Y(28)}
-2*gamma * Y(24);  % Y(29): hat{Y(29)}

];
    
end