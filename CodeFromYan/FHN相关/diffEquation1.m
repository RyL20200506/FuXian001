function f= diffEquation1(N)
    %% 生成微分方程
    %% input：
    %   N:节点数
    %   W：连接矩阵W[][]
    %% output：
    %   f
    f='@(v)[';
    for i=1:N
        if i==1
        str=['v(' num2str(i) ')-1/3* power(v(' num2str(i) '),3)-v(' num2str(N+i) ')+c1'];
        for j=1:N
            str=[str,'+W(' num2str(i) ',' num2str(j) ')*v(',num2str(j),')'];
        end

        else
        str=['v(' num2str(i) ')-1/3* power(v(' num2str(i) '),3)-v(' num2str(N+i) ')'];
        for j=1:N
            str=[str,'+W(' num2str(i) ',' num2str(j) ')*v(',num2str(j),')'];
        end
        end
        f=[f,str,';'];      
    end
    for i=1:N
         str1=['delta*(v(' num2str(i) ')-b*v(' num2str(N+i) ')+c2)']; 
         f=[f,str1,';'];
    end
    f=[f,']'];
end

