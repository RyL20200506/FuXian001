function [BASE] = baseEuqation(B,N,s)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%����һ�����󣬴������������Ϊ��λ��
% len=length(A);
% N��ʾ�м���������
% ����
% syms a1 a2 a3 a4
% B=[a1,a2,a3,a4];
% N=2;s=1;
[a,b]=size(B); %a��ʾÿ����������Ŀ��b��ʾ�м�������
    BASE=[];
    for i=1:a
        BASE1=[1;B(i,s);B(i,s+N)];
%         for j=s:N:b
%             BASE1=[BASE1;B(i,j)];
%         end       
        %����
        for k=s:N:b
             for p=k:N:b
                BASE1=[BASE1;B(i,k)*B(i,p)];
             end
        end
    %����
        for l=s:N:b
            for o=l:N:b
                for q=o:N:b
                    BASE1=[BASE1;B(i,l)*B(i,o)*B(i,q)];
                end
             end
        end
      %�Ĵ�
%         for l=s:N:b
%             for o=l:N:b
%                 for q=o:N:b
%                     for r=q:N:b
%                        BASE1=[BASE1;B(i,l)*B(i,o)*B(i,q)*B(i,r)];
%                     end
%                 end
%              end
%         end
%         for ss=1:N  
%             BASE1=[BASE1;B(i,ss)];
%         end
    BASE=[BASE,BASE1];
    end

   