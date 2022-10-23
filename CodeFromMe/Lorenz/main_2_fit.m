% 用于学习多元函数拟合

% 构造训练集
rng default % for reproducibility
tdata = 0:0.1:10;
ydata = 40*exp(-0.5*tdata) + randn(size(tdata));

% 构造需要最小化的函数, x是希望拟合的参数, 输入tdata进去, 希望得到ydata
fun = @(x)fun_2_sseval(x,tdata,ydata);

% 搜索最优参数
x0 = rand(2,1);
bestx = fminsearch(fun,x0);

%检查拟合质量
A = bestx(1);
lambda = bestx(2);
yfit = A*exp(-lambda*tdata);
plot(tdata,ydata,'*');
hold on
plot(tdata,yfit,'r');
xlabel('tdata')
ylabel('Response Data and Curve')
title('Data and Best Fitting Exponential Curve')
legend('Data','Fitted Curve')
hold off
