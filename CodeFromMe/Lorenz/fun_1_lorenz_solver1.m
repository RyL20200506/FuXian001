% �����Ȼ��� �������Ķ���ѧ, �ֶ����ΰ�
function ret = fun_1_lorenz_solver1(t,Y)

% x:y(1), y:y(2), z:y(3), 
% hat{x}: y(4), hat{y}: y(5), hat{z}: y(6)
% hat{a}: y(7), hat{b}: y(8), hat{r}: y(9), 
% e_1:y(10), e_2:y(11)
% G_k_a: y(12), G_k_b: y(13), G_k_r: y(14)
% G_k_e_1: y(15), G_k_e_2: y(16)
% H_k_a: y(17), H_k_b: y(18), H_k_r: y(19)
% H_k_e_1: y(20), H_k_e_2: y(21)
% ��ʽ�ο�:http://127.0.0.1:1111/center/call_exp_general_function/ftn.url_open_view_via_id('20220912-1853-4398-4377-000000000000')

t;  % ������֤ʱ����øú���ʱ�����õĲ���, ���ֲ����ǵȲ�����
% a=15;
% b=3;
% r=35;

a=10;
b=8/3;
r=28;

gamma=0.001;
alpha=2;
beta=0.5;
% ��дһ�ֱȽ�ͨ�õ�����
ret=[

% ��������: ��ʵϵͳ
a*(Y(2)-Y(1));  % Y(1):Y(1), a*(Y(2)-Y(1))
r*Y(1)-Y(2)-Y(1)*Y(3);  % Y(2):Y(2), r*Y(1)-Y(2)-Y(1)*Y(3)
Y(1)*Y(2)-b*Y(3);  % Y(3):Y(3), Y(1)*Y(2)-b*Y(3)

% ��������: ���ϵͳ
Y(27)*(Y(5)-Y(1));  % Y(4):Y(4), Y(27)*(Y(5)-Y(1))
Y(29)*Y(1)-Y(5)-Y(1)*Y(6)+Y(30);  % Y(5):Y(5)
Y(1)*Y(5)-Y(28)*Y(6)+Y(31);  % Y(6):Y(6)

% ��������: Y(1)ƫ������ض���ѧ
(Y(5)-Y(1))+Y(27)*Y(12);  % Y(7)
Y(27)*Y(13);  % Y(8)
Y(27)*Y(14);  % Y(9)
Y(27)*Y(15);  % Y(10)
Y(27)*Y(16);  % Y(11)

% ��������: Y(2)ƫ������ض���ѧ
-Y(12)+(-Y(1))*Y(17);  % Y(12):Y(12) = -Y(5)_a
-Y(13)+(-Y(1))*Y(18);  % Y(13):Y(13)
Y(1)-Y(14)+(-Y(1))*Y(19);  % Y(14):Y(14)
-Y(15)+(-Y(1))*Y(20)+1;  % Y(5)_hat{Y(30)}:Y(5)_hat{Y(30)}
-Y(16)+(-Y(1))*Y(21);  % Y(5)_hat{Y(31)}:Y(5)_hat{Y(31)}

% ��������: Y(3)ƫ������ض���ѧ
Y(1)*Y(12)-Y(28)*Y(17);  % Y(17):Y(17)
Y(1)*Y(13)+(-1)*(Y(6)+Y(28)*Y(18));  % Y(18):Y(18)
Y(1)*Y(14)-Y(28)*Y(19);  % Y(19):Y(19)
Y(1)*Y(15)-Y(28)*Y(20);  % Y(20)
Y(1)*Y(16)-Y(28)*Y(21)+1;  % Y(21)

% Deltaƫ���Ķ���ѧ
-alpha*Y(22) + (-2)*( (a*(Y(2)-Y(1)))-(Y(27)*(Y(5)-Y(1))) ) * ((Y(5)-Y(1)) + Y(27)*Y(12));  % Y(22)
-alpha*Y(23) + (-2)*( (a*(Y(2)-Y(1)))-(Y(27)*(Y(5)-Y(1))) ) * (Y(27)*Y(13));  % Y(23)
-alpha*Y(24) + (-2)*( (a*(Y(2)-Y(1)))-(Y(27)*(Y(5)-Y(1))) ) * (Y(27)*Y(14));   % Y(24)
% Delta����ei�Ķ���ѧ: �ܹؼ�!
-alpha*Y(25) + (-2)*( (a*(Y(2)-Y(1))) - (Y(27)*(Y(5)-Y(1))) ) * Y(27)*Y(15) + 2*(0.5)* (-2*gamma * Y(25));  % Y(25)  % Ӧ����д����.
-16.37*Y(26) + (-2)*( (a*(Y(2)-Y(1))) - (Y(27)*(Y(5)-Y(1))) ) * Y(27)*Y(16) + 2*(0.5)* (-2*gamma * Y(26));  % Y(26)

% ��������ѧ Ӧ������ 1
-2*gamma * Y(22)  % Y(27): Y(27)
-2*gamma * Y(23)  % Y(28): Y(28)
-2*gamma * Y(24)  % Y(29): Y(29)
-2*gamma * Y(25)  % Y(30): hat{Y(30)}
-2*gamma * Y(26)  % Y(31): hat{Y(31)}

];
    
end