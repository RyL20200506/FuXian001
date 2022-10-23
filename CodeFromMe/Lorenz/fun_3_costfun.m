% 复现: Wang - 重构动力学
% Wang 误差函数

function cost = fun_3_costfun(PARA, X, Y)
% PARA是希望拟合的参数
% X是特征
% Y是希望进行的预测

% 测试输入
% PARA = ones(48,1);
% X = [1,1,1; 2,2,2];
% Y = [3,3,3; 4,4,4];

% 正文
size_X = size(X);  % 获取训练特征的规格
pred_Y = [];  % 容器: 基于X得到的Y的预测

% 正文: 对各个X分别预测下一时刻的数值
for i = 1:size_X(1)  % 对训练特征遍历
    a_X = X(i, :);  % 拿出一行来
    a_pred_Y = fun_3_4_predict(PARA, a_X, 1);  % 1表示只预测未来1步
    pred_Y = [pred_Y; a_pred_Y];  % 新增一行
end

% 正文: 计算误差, 并输出
cost = sum( sum( (pred_Y-Y).^2 ) )   % 做差, 点乘平方, 对行求和, 对列求和





