function [BASE] = fun_16_baseEuqation1_with_para_4(B, para)
% para: 表示我们需要考虑几次

[a,b]=size(B); %a表示每个分量的数目，b表示有几个分量
n = para;  % 次数, 

% 算法: 

BASE=[];
for i=1:a   
    BASE1=[1];
    % %  一次
    for j=1:b
        BASE1=[BASE1;B(i,j)];
    end
% %     %二次
    for k=1:b
        for p=k:b
            BASE1=[BASE1;B(i,k)*B(i,p)];
        end
    end
% %     %三次
    for k=1:b
        for p=k:b
            for l=p:b
                BASE1=[BASE1;B(i,k)*B(i,p)*B(i,l)];
            end
        end
    end
% %     %四次
    for k=1:b
        for p=k:b
            for l=p:b
                for o=l:b
                    BASE1=[BASE1;B(i,k)*B(i,p)*B(i,l)*B(i,o)];
                end
            end
        end
    end
BASE=[BASE,BASE1];
end