function [BASE] = fun_16_baseEuqation1_with_para_4(B, para)
% para: ��ʾ������Ҫ���Ǽ���

[a,b]=size(B); %a��ʾÿ����������Ŀ��b��ʾ�м�������
n = para;  % ����, 

% �㷨: 

BASE=[];
for i=1:a   
    BASE1=[1];
    % %  һ��
    for j=1:b
        BASE1=[BASE1;B(i,j)];
    end
% %     %����
    for k=1:b
        for p=k:b
            BASE1=[BASE1;B(i,k)*B(i,p)];
        end
    end
% %     %����
    for k=1:b
        for p=k:b
            for l=p:b
                BASE1=[BASE1;B(i,k)*B(i,p)*B(i,l)];
            end
        end
    end
% %     %�Ĵ�
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