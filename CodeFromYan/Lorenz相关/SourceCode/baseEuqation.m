function [BASE] = baseEuqation(B)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%����һ�����󣬴������������Ϊ��λ��
% len=length(A);
[a,b]=size(B); %a��ʾÿ����������Ŀ��b��ʾ�м�������
BASE=[];
for i=1:a
    BASE1=[1];
    for j=1:b
        BASE1=[BASE1;B(i,j)];

    end
    %����
    for k=1:b
        for p=k:b
            BASE1=[BASE1;B(i,k)*B(i,p)];
        end
    end
    %����
    for l=1:b
        for o=l:b
            for q=o:b
                BASE1=[BASE1;B(i,l)*B(i,o)*B(i,q)];
            end
        end
    end
BASE=[BASE,BASE1];
end


