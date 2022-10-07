function [BASE] = baseEuqation(B,N,s)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%传入一个矩阵，传出以这个矩阵为单位的
% len=length(A);
% N表示有几个网络格点
% 测试
% syms a1 a2 a3 a4
% B=[a1,a2,a3,a4];
% N=2;s=1;
[a,b]=size(B); %a表示每个分量的数目，b表示有几个分量
    BASE=[];
    for i=1:a
        BASE1=[1;B(i,s);B(i,s+N)];
%         for j=s:N:b
%             BASE1=[BASE1;B(i,j)];
%         end       
        %二次
        for k=s:N:b
             for p=k:N:b
                BASE1=[BASE1;B(i,k)*B(i,p)];
             end
        end
    %三次
        for l=s:N:b
            for o=l:N:b
                for q=o:N:b
                    BASE1=[BASE1;B(i,l)*B(i,o)*B(i,q)];
                end
             end
        end
      %四次
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

   