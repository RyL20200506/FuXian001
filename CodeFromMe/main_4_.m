% 线性拟合的demo
% 参考: https://ww2.mathworks.cn/help/stats/regress.html

% 首先构建基向量, 然后线性拟合求解
clc; clear

load carsmall
x1 = Weight;
x2 = Horsepower;    % Contains NaN data
y = MPG;  % 这个就是希望得到的预测结果.

% 构建向量
X = [ones(size(x1)) x1 x2 x1.*x2];  
b = regress(y,X)    % Removes NaN data % 其实还会有其他参数, 比如置信度和残差之类的


% 绘制散点图
scatter3(x1,x2,y,'filled')
hold on
x1fit = min(x1):100:max(x1);
x2fit = min(x2):10:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('Weight')
ylabel('Horsepower')
zlabel('MPG')
view(50,10)
hold off
