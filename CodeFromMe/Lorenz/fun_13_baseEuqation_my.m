
% ��֮ǰ������, ��߾�������Ч��
function [BASE1] = fun_13_baseEuqation_my(B)

% test_input
% B = [2, 3, 4; 1,2,3];

% Content:
[a,b]=size(B); % a��ʾÿ����������Ŀ��b��ʾ�м�������

BASE1 = ones(a,1);
for j=1:b 
    BASE1=[BASE1, B(:,j)];  % ����ϵ��x, y, z
end
%����
for k=1:b
    for p=k:b
        BASE1=[BASE1, B(:,k).*B(:,p)];  % ���Ƕ����� xx xy xz yy yz zz
    end
end
%����
for l=1:b
    for o=l:b
        for q=o:b
            BASE1=[BASE1, B(:,l).*B(:,o).*B(:,q)];  % xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
        end
    end
end

end


