function [BASE] = fun_16_baseEuqation1_with_para_10(B)
% para: 表示我们需要考虑几次

[a,b]=size(B); %a表示每个分量的数目，b表示有几个分量
% n = para;  % 次数, 

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
%     %五次
    for k=1:b
        for p=k:b
            for l=p:b
                for o=l:b
                    for m=o:b
                        BASE1=[BASE1;B(i,k)*B(i,p)*B(i,l)*B(i,o)*B(i,m)];
                    end
                end
            end
        end
    end
% %     %6次
%     for k1=1:b
%         for k2=k1:b
%             for k3=k2:b
%                 for k4=k3:b
%                     for k5=k4:b
%                         for k6=k5:b
%                             BASE1=[BASE1;B(i,k1)*B(i,k2)*B(i,k3)*B(i,k4)*B(i,k5)*B(i,k6)];
%                         end
%                     end
%                 end
%             end
%         end
%     end    
% %     %7次
%     for k1=1:b
%         for k2=k1:b
%             for k3=k2:b
%                 for k4=k3:b
%                     for k5=k4:b
%                         for k6=k5:b
%                             for k7=k6:b
%                                 BASE1=[BASE1;B(i,k1)*B(i,k2)*B(i,k3)*B(i,k4)*B(i,k5)*B(i,k6)*B(i,k7)];
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end    
% %     %8次
%     for k1=1:b
%         for k2=k1:b
%             for k3=k2:b
%                 for k4=k3:b
%                     for k5=k4:b
%                         for k6=k5:b
%                             for k7=k6:b
%                                 for k8=k7:b
%                                     BASE1=[BASE1;B(i,k1)*B(i,k2)*B(i,k3)*B(i,k4)*B(i,k5)*B(i,k6)*B(i,k7)*B(i,k8)];
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end    
BASE=[BASE,BASE1];
end