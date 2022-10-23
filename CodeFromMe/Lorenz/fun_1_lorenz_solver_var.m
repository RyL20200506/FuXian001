% 洛伦兹混沌 求解参数的动力学
% 区别: 通过变量名来完成的
function ret = fun_1_lorenz_solver(t,Y)

t;  % 用于验证时间调用该函数时所是用的步长, 发现并不是等步长的

a=10;
b=8/3;
r=28;
% a=15;
% b=3;
% r=35;

gamma=0.001;
alpha=3;
beta=0.5;

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

% 下面的还没改

% 先写一种比较通用的语言
ret=[

% 绝对无误: 真实系统
a*(y-x);  % x:x, a*(y-x)
r*x-y-x*z;  % y:y, r*x-y-x*z
x*y-b*z;  % z:z, x*y-b*z
% 绝对无误: 拟合系统
hata * (haty-x);  % hat{x}:hat{x}, hat{a}*(hat{y}-x)
hatr*x - haty - x*hatz+e1;  % hat{y}:hat{y}
x*haty - hatb*hatz + e2;  % hat{z}:hat{z}
% 绝对无误: x偏导的相关动力学
(hat{y}-x)+hat{a}*hat{y}hat{a};  % hat{x}_hat{a}
hat{a}*hat{y}_hat{b};  % hat{x}_hat{b}
hat{a}*hat{y}_hat{r};  % hat{x}_hat{r}
hat{a}*hat{y}_e1;  % hat{x}_e1
hat{a}*hat{y}_e2;  % hat{x}_e2
% 绝对无误: y偏导的相关动力学
-hat{y}_hat{a}+(-x)*hat{z}_hat{a};  % hat{y}_hat{a}:hat{y}_hat{a} = -hat{y}_a
-hat{y}_hat{b}+(-x)*hat{z}_hat{b};  % hat{y}_hat{b}:hat{y}_hat{b}
x-hat{y}_hat{r}+(-x)*hat{z}_hat{r};  % hat{y}_hat{r}:hat{y}_hat{r}
-hat{y}_e1+(-x)*hat{z}_e1+1;  % hat{y}_hat{e1}:hat{y}_hat{e1}
-hat{y}_e2+(-x)*hat{z}_e2;  % hat{y}_hat{e2}:hat{y}_hat{e2}
% 绝对无误: z偏导的相关动力学
x*hat{y}_hat{a}-hat{b}*hat{z}_hat{a};  % hat{z}_hat{a}:hat{z}_hat{a}
x*hat{y}_hat{b}+(-1)*(hat{z}+hat{b}*hat{z}_hat{b});  % hat{z}_hat{b}:hat{z}_hat{b}
x*hat{y}_hat{r}-hat{b}*hat{z}_hat{r};  % hat{z}_hat{r}:hat{z}_hat{r}
x*hat{y}_e1-hat{b}*hat{z}_e1;  % hat{z}_e1
x*hat{y}_e2-hat{b}*hat{z}_e2+1;  % hat{z}_e2
% Delta偏导的动力学
-alpha*G_hat{a} + (-2)( (a*(y-x))-(hat{a}*(hat{y}-x)) ) * ((hat{y}-x) + hat{a}*hat{y}_hat{a});  % G_hat{a}
-alpha*G_hat{b} + (-2)( (a*(y-x))-(hat{a}*(hat{y}-x)) ) * (hat{a}*hat{y}_hat{b});  % G_hat{b}
-alpha*G_hat{r} + (-2)( (a*(y-x))-(hat{a}*(hat{y}-x)) ) * (hat{a}*hat{y}_hat{r});   % G_hat{r}
% 下面的最重要
-alpha*G_e1 + (-2)( (a*(y-x)) - (hat{a}*(hat{y}-x)) ) * hat{a}*hat{y}_e1 + 2*beta* (-2*gamma * G_e1);  % G_e1
-alpha*G_e2 + (-2)( (a*(y-x)) - (hat{a}*(hat{y}-x)) ) * hat{a}*hat{y}_e2 + 2*beta* (-2*gamma * G_e2);  % G_e2
% 参数动力学 应该无误 1
-2*gamma * G_hat{a}  % hat{a}: hat{a}
-2*gamma * G_hat{b}  % hat{b}: hat{b}
-2*gamma * G_hat{r}  % hat{r}: hat{r}
-2*gamma * G_e1  % e1: hat{e1}
-2*gamma * G_e2  % e2: hat{e2}

];
    
end