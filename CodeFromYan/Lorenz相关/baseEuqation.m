function [BASE] = baseEuqation(B)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%传入一个矩阵，传出以这个矩阵为单位的
% len=length(A);

% test
% B = [[1,2,3]; [4,5,6]]  % 样例: 用于理解

[a,b] = size(B);  % R: a: 行数|分量维数，b: 列数|表示分量个数
BASE=[];
for i=1:a  % 对每一行, 比如第1行, 即第1个时刻, 此时随着行的增加, 时间增加.
    BASE1=[1];  % 初始化常数项
    for j=1:b  % 对每一列, 对第1行的每一列: 
        BASE1=[BASE1;B(i,j)];  % 把元素合并过去, 注意这里是分号, 即放在这个元素的下面
    end
    
    %二次: R: 考虑第1个向量内部之间的两两相互作用
    for k=1:b
        for p=k:b  % 注意这里是从k开始往后了, 所以不会出现重复考虑
            BASE1=[BASE1;B(i,k)*B(i,p)] ;
        end
    end
    %三次: 同上, 是考虑第1个向量内部元素的两两相互作用.
    for l=1:b
        for o=l:b
            for q=o:b
                BASE1=[BASE1;B(i,l)*B(i,o)*B(i,q)];
            end
        end
    end

BASE=[BASE,BASE1];  % 对于每个列向量, 都进行同样的处理, 然后合并在一起. 进行输出.
end


