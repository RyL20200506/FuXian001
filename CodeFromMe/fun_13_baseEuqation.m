function [BASE] = fun_13_baseEuqation(B)

% test_input
% B = [2, 3, 4; 1,2,3];

% Content:
[a,b]=size(B); % a��ʾÿ����������Ŀ��b��ʾ�м�������

BASE=[];
for i=1:a  % ��B��ÿһ�����ɶ�Ӧ������
    BASE1=[1];
    for j=1:b 
        BASE1=[BASE1, B(i,j)];  % ����ϵ��x, y, z
    end
    %����
    for k=1:b
        for p=k:b
            BASE1=[BASE1, B(i,k)*B(i,p)];  % ���Ƕ����� xx xy xz yy yz zz
        end
    end
    %����
    for l=1:b
        for o=l:b
            for q=o:b
                BASE1=[BASE1, B(i,l)*B(i,o)*B(i,q)];  % xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
            end
        end
    end
    % ���кϲ�
    BASE=[BASE; BASE1];  
end


