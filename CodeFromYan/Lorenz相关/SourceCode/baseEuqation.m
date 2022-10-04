function [BASE] = baseEuqation(B)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%传入一个矩阵，传出以这个矩阵为单位的
% len=length(A);
[a,b]=size(B); %a表示每个分量的数目，b表示有几个分量
BASE=[];
for i=1:a
    BASE1=[1];
    for j=1:b
        BASE1=[BASE1;B(i,j)];

    end
    %二次
    for k=1:b
        for p=k:b
            BASE1=[BASE1;B(i,k)*B(i,p)];
        end
    end
    %三次
    for l=1:b
        for o=l:b
            for q=o:b
                BASE1=[BASE1;B(i,l)*B(i,o)*B(i,q)];
            end
        end
    end
BASE=[BASE,BASE1];
end


