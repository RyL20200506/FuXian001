% 基于差分法实现的动力学

% 拟合动力学: 需要外部输入差分
function ret = fun_6_5_lorenz_solver(t,Y)

% x:y(1), y:y(2), z:y(3), 
% hat{x}: y(4), hat{y}: y(5), hat{z}: y(6)
% hat{a}: y(7), hat{b}: y(8), hat{r}: y(9), 
% e_1:y(10), e_2:y(11)
% G_k_a: y(12), G_k_b: y(13), G_k_r: y(14)
% G_k_e_1: y(15), G_k_e_2: y(16)
% H_k_a: y(17), H_k_b: y(18), H_k_r: y(19)
% H_k_e_1: y(20), H_k_e_2: y(21)
% 公式参考:http://127.0.0.1:1111/center/call_exp_general_function/ftn.url_open_view_via_id('20220912-1853-4398-4377-000000000000')

% 准备:插值获取dot_x
load mat_6_dot_x.mat  % 导入离散数据
step_length = 0.00001;  % 步长
time_range = [0:step_length:20-step_length];  % 获取离散步长
% % vq1 = interp1(time_range,dot_x',[20,21]);
% F = griddedInterpolant(time_range,dot_x');

load fun_6_6_F

% 正确的代码
if (mod(floor(t), 10)==1)
    t
end
a=10;
b=8/3;
r=28;
gamma=0.002;
alpha=3;
beta=1;
% 先写一种比较通用的语言
ret=[

% 原系统
a*(Y(2)-Y(1));  % Y(1):Y(1), a*(Y(2)-Y(1))
r*Y(1)-Y(2)-Y(1)*Y(3);  % Y(2):Y(2), r*Y(1)-Y(2)-Y(1)*Y(3)
Y(1)*Y(2)-b*Y(3);  % Y(3):Y(3), Y(1)*Y(2)-b*Y(3)

% 拟合系统
Y(25)*(Y(5)-Y(1));  % Y(4):Y(4), Y(25)*(Y(5)-Y(1))
Y(27)*Y(1)-Y(5)-Y(1)*Y(6)+Y(28);  % Y(5):Y(5), 
Y(1)*Y(5)-Y(26)*Y(6)+Y(29);  % Y(6):Y(6)


1*(Y(5)-Y(1))+Y(25)*Y(10); % Y(7):Y(7)
Y(25)*Y(11);  % Y(8):Y(4)_b
Y(25)*Y(12);  % Y(9):Y(4)_r = Y(25) * hat(Y(2))_r

-Y(10)+(-Y(1))*Y(15);  % Y(10):Y(10) = -Y(5)_a
-Y(11)+(-Y(1))*Y(16);  % Y(11):Y(11)
Y(1)-Y(12)+(-Y(1))*Y(17);  % Y(12):Y(12)
-Y(13)+(-Y(1))*Y(18)+1;  % Y(5)_hat{Y(28)}:Y(5)_hat{Y(28)}
-Y(14)+(-Y(1))*Y(19);  % Y(5)_hat{Y(29)}:Y(5)_hat{Y(29)}

% 感觉无误
Y(1)*Y(10)-Y(26)*Y(15);  % Y(15):Y(15)
Y(1)*Y(11)+(-1)*(Y(6)+Y(26)*Y(16));  % Y(16):Y(16)
Y(1)*Y(12)-Y(26)*Y(17);  % Y(17):Y(17)
Y(1)*Y(13)-Y(26)*Y(18);  % Y(6)_hat{Y(28)}:Y(6)_hat{Y(28)}
Y(1)*Y(14)-Y(26)*Y(19)+1;  % Y(6)_hat{Y(29)}:Y(6)_hat{Y(29)}

% % 参数的动力学 感觉无误
% -alpha*Y(20)+(-2)*( interp1(time_range,dot_x',t) -(Y(25)*(Y(5)-Y(1))) )*(2*Y(25)*Y(10)+(Y(5)-Y(1)));
% -alpha*Y(21)+(-2)*( interp1(time_range,dot_x',t)-(Y(25)*(Y(5)-Y(1))) )*(2*Y(25)*Y(11));
% -alpha*Y(22)+(-2)*( interp1(time_range,dot_x',t)-(Y(25)*(Y(5)-Y(1))) )*(2*Y(25)*Y(12));
% % 误差梯度动力学 感觉无误
% -alpha*Y(23)+(-2)*( interp1(time_range,dot_x',t) - (Y(25)*(Y(5)-Y(1))) )*Y(25)*Y(13)+2*beta* (-2*gamma * Y(23));
% -alpha*Y(24)+(-2)*( interp1(time_range,dot_x',t) - (Y(25)*(Y(5)-Y(1))) )*Y(25)*Y(14)+2*beta* (-2*gamma * Y(24));

% 参数的动力学 感觉无误
-alpha*Y(20)+(-2)*( F(t) -(Y(25)*(Y(5)-Y(1))) )*(2*Y(25)*Y(10)+(Y(5)-Y(1)));
-alpha*Y(21)+(-2)*( F(t)-(Y(25)*(Y(5)-Y(1))) )*(2*Y(25)*Y(11));
-alpha*Y(22)+(-2)*( F(t)-(Y(25)*(Y(5)-Y(1))) )*(2*Y(25)*Y(12));
% 误差梯度动力学 感觉无误
-alpha*Y(23)+(-2)*( F(t) - (Y(25)*(Y(5)-Y(1))) )*Y(25)*Y(13)+2*beta* (-2*gamma * Y(23));
-alpha*Y(24)+(-2)*( F(t) - (Y(25)*(Y(5)-Y(1))) )*Y(25)*Y(14)+2*beta* (-2*gamma * Y(24));

% 参数动力学 应该无误
-2*gamma * Y(20);  % Y(25): Y(25)
-2*gamma * Y(21);  % Y(26): Y(26)
-2*gamma * Y(22);  % Y(27): Y(27)
-2*gamma * Y(23);  % Y(28): hat{Y(28)}
-2*gamma * Y(24);  % Y(29): hat{Y(29)}

];


% 查看误差的大小
r1 = a*(Y(2)-Y(1)) % 真实值
in1 = interp1(time_range, dot_x',t)  % 
f1 = F(t)  % F的插值结果
    
end