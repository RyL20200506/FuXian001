
% 在之前基础上, 提高矩阵运算效率
function [BASE1] = fun_13_baseEuqation_my(B)

% test_input
% B = [2, 3, 4; 1,2,3];

% Content:
[a,b]=size(B); % a表示每个分量的数目，b表示有几个分量

BASE1 = ones(a,1);
for j=1:b 
    BASE1=[BASE1, B(:,j)];  % 考虑系数x, y, z
end
%二次
for k=1:b
    for p=k:b
        BASE1=[BASE1, B(:,k).*B(:,p)];  % 考虑二阶项 xx xy xz yy yz zz
    end
end
%三次
for l=1:b
    for o=l:b
        for q=o:b
            BASE1=[BASE1, B(:,l).*B(:,o).*B(:,q)];  % xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
        end
    end
end

end


