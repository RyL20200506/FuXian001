function [BASE] = fun_13_baseEuqation(B)

% test_input
% B = [2, 3, 4; 1,2,3];

% Content:
[a,b]=size(B); % a表示每个分量的数目，b表示有几个分量

BASE=[];
for i=1:a  % 对B的每一行生成对应的数据
    BASE1=[1];
    for j=1:b 
        BASE1=[BASE1, B(i,j)];  % 考虑系数x, y, z
    end
    %二次
    for k=1:b
        for p=k:b
            BASE1=[BASE1, B(i,k)*B(i,p)];  % 考虑二阶项 xx xy xz yy yz zz
        end
    end
    %三次
    for l=1:b
        for o=l:b
            for q=o:b
                BASE1=[BASE1, B(i,l)*B(i,o)*B(i,q)];  % xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
            end
        end
    end
    % 将行合并
    BASE=[BASE; BASE1];  
end


