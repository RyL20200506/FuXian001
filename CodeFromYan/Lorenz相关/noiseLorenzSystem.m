% Taylor cacluator nosie Lorenz System
close all;clear;clc;  % R: 正常操作
% for D=1:3  


for zkkk=1:2  % ?猜: 应该是循环的次数 

    % C.1 参数声明
    delta=10;b=8/3;r=28;  % Lorenz正常参数
    N=10000;
    x0 = [-8;7;27];  % 初始值
    X_n = [x0];  % R: 容器
    t0 = 0.01;  % 单位计算时间
    D = 500;  % 噪声方差
    m = 1;  % R: 哦这个好像是时间窗口之类的

    % C.2. 计算微分方程的轨迹 并加上噪声
    b1=sqrt(2*D*t0); % ??? 这是啥, 反正跟噪声有关
    b2=1/2*b1*t0;
    b3=1/sqrt(3)*b2;   
    w1=randn(N,3);w2=randn(N,3);  % R: 生成标准白噪声
    for i=1:N  
        s1=b1*w1(i,:);  % 
        s2=b2*w1(i,:)+b3*w2(i,:);  % 
        df2=[0;(r-X_n(3,i)-1-X_n(1,i));(X_n(2,i)+X_n(1,i)-b)];  %

        % R: 下面是Lorenz方程, 应该是原方程, 而不是拟合方程.
        fx=delta*(X_n(2,i)-X_n(1,i));  %
        fy=r*X_n(1,i)-X_n(2,i)-X_n(1,i)*X_n(3,i);  %
        fz=X_n(1,i)*X_n(2,i)-b*X_n(3,i);

        % R: 总之应该是往后推一段X_n
        df1=[-delta*fx+delta*fy; (r-X_n(3,i))*fx-fy-X_n(1,i)*fz; X_n(2,i)*fx+X_n(1,i)*fy-b*fz];  % 
        X_t=X_n(:,i)+[fx;fy;fz]*t0+1/2*t0*t0*df1+df2.*s2'+s1';  % R: 用最后一列数据去更新状态, 而且在状态更新时引入了噪声.
        X_n=[X_n,X_t];  % 合并新的列, 最后得到一个3行, 10001列的矩阵. 随着列的增加时间增加.
    end

    % C.3 画出计算得到的轨迹图
    A=X_n(:,1:end-1);  % 去掉了最后一列, 因为最后一列只作为结果判断, 不作为已知
    subplot(1,2,1)  % 
    plot(A(1,:),A(3,:),'linewidth',1);
    xlabel('x');
    ylabel('z');
    set(gca,'FontSize',30) ;
    subplot(1,2,2)
    plot(A(1,:),'linewidth',1);
    xlabel('step');
    ylabel('x');
    set(gca,'FontSize',30) ;

    % C.4 获取拟合所使用的基向量
    BASE = [];  % ? 感觉这句话可以删掉
    base = baseEuqation(A');  % A是随着列增加时间增加, A'修改为随着行增加时间增加.
    BASE = [base*t0;base*t0*t0];  % R: 这里t0对应原文的t_{unit}的概念, 其实就是g_i和h_i, 它们的基向量相同, 它们的
    
    % C.5 (Yan)m=1,从第一个元素算起,初始点数为1：L，对2：L+1拟合的结果
    DBASE=[base;2*t0*base];  % R 以另一种形式构造了基向量, 用的是原基向量和一阶基向量
    BASE=BASE(:,1:N-1);  % 去掉最后一列基向量 ???但是为什么要这样做, 因为BASE是基于base, base是基于A, A已经从X_n中去掉最后一行列.
    BASE_pinv=pinv(BASE);  % 准备好基向量的伪逆

    % ??? 为什么: 是根据: 差分=系数*BASE, 所以系数=差分/BASE吗?, 那diff_x不就是系数的含义吗
    diff_x=(A(1,2:end)-A(1,1:end-1))*BASE_pinv;  % R 做差, 再除以基向量, 这样仿佛能得到, 一个向量使得:这个向量*基向量=与现在的差距.
    diff_y=(A(2,2:end)-A(2,1:end-1))*BASE_pinv;  % 
    diff_z=(A(3,2:end)-A(3,1:end-1))*BASE_pinv;  % 
    
    % ??? DBASE的作用是什么, DX表示的是导数的含义吗? 为什么可以这样来计算导数?
    DX=diff_x*DBASE;  % 疑:但BASE和DBASE是两种不同的基向量.
    DY=diff_y*DBASE;  % 但这个变量之后才会用到, 我也不知道有没有用
    DZ=diff_z*DBASE;

    % ??? 这难道是得到了下一时刻的状态矩阵, 所以diff_x确实表示系数的含义, 因为diff_x*BASE就得到了差分.
    % ??? 但这也没有预测呀, 因为A(1,2:end)不是属于已知的部分吗?
    X=A(1,1:end-1)+diff_x*BASE;  % R 得到拟合后x分量的数据
    Y=A(2,1:end-1)+diff_y*BASE;  % 
    Z=A(3,1:end-1)+diff_z*BASE;  % 

    % X=diff_x*BASE
    % Y=diff_y*BASE;
    % Z=diff_z*BASE;

    % plot(X_n(1,:),X_n(2,:))
    
    
    for i=1:N  % ??? A应该只有3行, 下面就是第2行减第1行, 相当于y-x, 计算每一个时刻y-x的值.
        base_a(i,:)=(A(2,i)-A(1,i));  % ??? 这不是y-x吗, 跟a有什么关系呢? base_a怎么来理解?
    end
    for i=1:N
        base_b(i,:)=A(3,i);  % ??? 完全照抄z
    end
    for i=1:N
        base_r(i,:)=A(1,i);  % ??? 完全照抄x, 但为什么叫_r, x与r有什么关系吗?
    end
%     for i=1:9999
%         base_a(i,:)=(Y(1,i)-X(1,i))*t0;
%     end
%     for i=1:9999
%         base_b(i,:)=Z(1,i)*t0;
%     end
%     for i=1:9999
%         base_r(i,:)=X(1,i)*t0;
%     end
    
    base_ap=pinv(base_a);  % ??? 对y-x求逆?
    base_bp=pinv(base_b);  % 
    base_rp=pinv(base_r);
    
    for i=1:N
        cha_x(i,:)=DX(1,i);  % Why?, 为什么突然用DX去了
    end
    for i=1:N
        cha_z(i,:)=A(1,i)*A(2,i)-DZ(1,i);  % 为什么要这样写
    end
      for i=1:N
            cha_y(i,:)=DY(1,i)+A(2,i)+A(1,i)*A(3,i);  % 不懂啊不懂, 这代码真的可以用吗
      end
      a_value(zkkk,1)=base_ap*cha_x;  % a_value记录就是拟合的a的历史.
      b_value(zkkk,1)=base_bp*cha_z;  % b_value
      r_value(zkkk,1)=base_rp*cha_y;
end

a_value(end)
b_value(end)
r_value(end)

stop
%     %用论文11式精确方程
%     for i=1:size(A,2)
%         base_a(i,:)=(A(2,i)-A(1,i))*t0;  % C: base_a就是对A的两个相邻的时间做差所得;
%     end
%     base_ap=pinv(base_a); 
%     BASE_11=[base*t0; base*t0*t0];  % BASE_11跟BASE很像, 只是BASE少了一个维度而已;
% 
%     % a_value记录预测值, 问: 但预测的的完整公式是什么, 怎么跟论文中的理论对应上?
%     for i=1:size(A,2)
%        DERA(i,:)=0.5*t0*t0*(-a_value(zkkk,1)*a_value(zkkk,1)*(A(2,i)-A(1,i))+a_value(zkkk,1)*(r_value(zkkk,1)*A(1,i)-A(2,i)-A(1,i)*A(3,i)));
%     end
%     a_value(zkkk,1)=base_ap*(BASE_11'*diff_x'-DERA);  % C: 这里记录新的a值
% 
%     % 更新r
%     for i=1:size(A,2)
%         base_y(i,:)=A(1,i)*t0;  % 为啥?
%     end
%     base_yp=pinv(base_y);
%     for i=1:size(A,2)
%        DERA_y(i,:)=0.5*t0*t0*( a_value(zkkk,1)*(A(2,i)-A(1,i))*(r_value(zkkk,1)-A(3,i))-( r_value(zkkk,1)*A(1,i)-A(2,i)-A(1,i)*A(3,i))-A(1,i)*(A(2,i)*A(1,i)-b_value(zkkk,1)*A(3,i)));
%     end
%     r_value(zkkk,1)=base_yp*(BASE_11'*diff_y'-DERA_y+t0*(A(2,:)+A(1,:).*A(3,:))');  % C: 这里记录新的r值
% 
%     % 更新b
%     for i=1:size(A,2)
%         base_z(i,:)=A(3,i)*t0;
%     end
%     base_zp=pinv(base_z);
%     for i=1:size(A,2)
%        DERA_z(i,:)=0.5*t0*t0*(a_value(zkkk,1)*(A(2,i)-A(1,i))*A(2,i)+(r_value(zkkk,1)*A(1,i)-A(2,i)-A(1,i)*A(3,i))*A(1,i)-b_value(zkkk,1)*(A(2,i)*A(1,i)-b_value(zkkk,1)*A(3,i)));
%     end
%     b_value(zkkk,1)=-base_zp*(BASE_11'*diff_z'-DERA_z-t0*(A(1,:).*A(2,:))');
% 
%     % 计算误差
%     err_a(zkkk,1)=abs(a_value(zkkk,1)-delta)/delta;  % a也有另外一个名字叫delta
%     err_b(zkkk,1)=abs(b_value(zkkk,1)-b)/b; 
%     err_r(zkkk,1)=abs(r_value(zkkk,1)-r)/r;  % 用最新的(位于第zkkk行的)r_value与真实的r去做相对误差
%     noise_DD=0;
% 
%     for i=1:N-1
%        noise_D=(X(1,i)-A(1,i+1))^2;
%        noise_DD= noise_DD+ noise_D;
%     end
% 
%     noise_DD= noise_DD/(N-1);  % C 
%     D_D(zkkk,1)=noise_DD/(2*t0);  % C 对噪声的评估
%     err_D(zkkk,1)=abs(D_D(zkkk,1)-D)/D;  % C 
    % delta=a_value;b=b_value;r=-r_value;
% end
% 
% a_err=(mean(a_value)-delta)/delta;
% b_err=(mean(b_value)-b)/b;
% r_err=(mean(r_value)-r)/r;
% D_err=(mean(D_D)-D)/D;
% 
% a_err=mean(err_a)  % C: 对所有误差取均值, 这应该是最后的结果了
% b_err=mean(err_b)  
% r_err=mean(err_r)
% D_err=mean(err_D)  

%  N=10000;
% % x0=[-1;3;4];
% X_n=[];
% x0=[-8;7;27];
% 
% X_n=[x0];
% t0=0.01;%单位计算时间
% % D=1;%噪声方差
% D=0;
% b1=sqrt(2*D*t0);b2=1/2*b1*t0;b3=1/sqrt(3)*b2;
% w1=randn(N,3);w2=randn(N,3);
% for i=1:N
%     s1=b1*w1(i,:);
%     s2=b2*w1(i,:)+b3*w2(i,:);
%     df2=[0;(r-X_n(3,i)-1-X_n(1,i));(X_n(2,i)+X_n(1,i)-b)];
%     fx=delta*(X_n(2,i)-X_n(1,i));
%     fy=r*X_n(1,i)-X_n(2,i)-X_n(1,i)*X_n(3,i);
%     fz=X_n(1,i)*X_n(2,i)-b*X_n(3,i);
%     df1=[-delta*fx+delta*fy;(r-X_n(3,i))*fx-fy-X_n(1,i)*fz;X_n(2,i)*fx+X_n(1,i)*fy-b*fz]; 
%     X_t=X_n(:,i)+[fx;fy;fz]*t0+1/2*t0*t0*df1+df2.*s2'+s1';
%     X_n=[X_n,X_t];
% end
% AA=X_n(:,1:N);
% end
% hold on
% plot(X,Y)
