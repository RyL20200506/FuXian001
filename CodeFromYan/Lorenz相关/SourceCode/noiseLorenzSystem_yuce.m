%% Taylor cacluator nosie Lorenz System
close all;clear;clc;
% 
delta=10;b=8/3;r=28;
N=15000;
% x0=[-1;3;4];
x0=[-8;7;27];
X_n=[x0];
t0=0.01;%单位计算时间
m=1;
D=0;%噪声方差
% D=0;
b1=sqrt(2*D*t0);b2=1/2*b1*t0;b3=1/sqrt(3)*b2;
w1=randn(N,3);w2=randn(N,3);
for i=1:N
    s1=b1*w1(i,:);
    s2=b2*w1(i,:)+b3*w2(i,:);
    df2=[0;(r-X_n(3,i)-1-X_n(1,i));(X_n(2,i)+X_n(1,i)-b)];
    fx=delta*(X_n(2,i)-X_n(1,i));
    fy=r*X_n(1,i)-X_n(2,i)-X_n(1,i)*X_n(3,i);
    fz=X_n(1,i)*X_n(2,i)-b*X_n(3,i);
    df1=[-delta*fx+delta*fy;(r-X_n(3,i))*fx-fy-X_n(1,i)*fz;X_n(2,i)*fx+X_n(1,i)*fy-b*fz]; 
    X_t=X_n(:,i)+[fx;fy;fz]*t0+1/2*t0*t0*df1+df2.*s2'+s1';
    X_n=[X_n,X_t];
end

% % Rossler 系统
% x0=[4;-4;5];
% a=0.2;b=0.2;c=5.7;
% N=4000;
% X_n=[x0];
% t0=0.01;%单位计算时间
% % D=0.5;%噪声方差
% m=1;
% D=0;
% b1=sqrt(2*D*t0);b2=1/2*b1*t0;b3=1/sqrt(3)*b2;
% w1=randn(N,3);w2=randn(N,3);
% for i=1:N
%     s1=b1*w1(i,:);
%     s2=b2*w1(i,:)+b3*w2(i,:);
%     df2=[-2;1+a;X_n(3,i)+X_n(1,i)-c];
%     fx=-X_n(2,i)-X_n(3,i);
%     fy=X_n(1,i)+a*X_n(2,i);
%     fz=b+X_n(3,i)*(X_n(1,i)-c);
%     df1=[-1*fy-1*fz;fx+a*fy;X_n(3,i)*fx+(X_n(1,i)-c)*fz]; 
%     X_t=X_n(:,i)+[fx;fy;fz]*t0+1/2*t0*t0*df1+df2.*s2'+s1';
%     X_n=[X_n,X_t];
% end




B=X_n(:,1:12000);

k_length=length(B);
%%%%重构维度的构建
% BX=abs(BX);
% M_ZHI=1;
M=3;
% tau=16;
% % t0=t;
% A=reconstitution(B,k_length,M,tau);
A=B;
% J=1;
% for i=1:10000
%     A(2,J)=A(1,J*3);
%     A(3,J)=A(1,J*6);
%     J=J+1;
% end
% a=5;
%  A(1,:)=B(3,13:10000);
%  A(2,:)=B(3,7:9994);
%  A(3,:)=B(3,1:9988);

%  A(1,:)=B(1,21:10000);
%  A(2,:)=B(1,11:9990);
%  A(3,:)=B(1,1:9980);
%  
%%重构维度的构建
% A=[];
% q=1;
% tau=2;
% M=3;
% for i=1:M
%     A=[A,B(2,1+tau*(q-1):size(B,2)-tau*(M-q))'];
%     q=q+1;
% end
% A=A';
%  A(1,:)=B(1,1:9980);
%  A(2,:)=B(1,11:9990);
%  A(3,:)=B(1,21:10000);
%  A(1,:)=B(1,1:9990);
%  A(2,:)=B(1,11:10000);
%  A(3,:)=B(3,1:9990);

%  A(1,:)=B(3,1:9988);
%  A(2,:)=B(3,7:9994);
%  A(3,:)=B(3,13:10000);


% A=A(:,1:10000);
%对A矩阵svd分解，找前三个最大特征值对应的特征向量
% [UU,SS,VV] = svd(A);
% A=VV(:,1:M)';
% VV=VV';
% figure(1);
% plot(VV(:,1),VV(:,2))
% figure(2);
% plot(VV(:,1),VV(:,3))
% figure(3);
% plot(VV(:,2),VV(:,3))
% figure(4);
% plot3(VV(:,1),VV(:,2),VV(:,3))
BASE=[];
DBASE=[];
A_VAL=[];
XA_VAL=[];
YA_VAL=[];
ZA_VAL=[];
A_index=1;
for i=m+1:size(A,2)-m
    for j=-m:m
       if j~=0
     
       base=baseEuqation(A(:,i)')';
       BASE=[BASE;base*j*t0,base*j*t0*j*t0];%m=1,从第一个元素算起,初始点数为1：L，对2：L+1拟合的结果
       DBASE=[DBASE;base,2*base*j*t0];
%        for i_m=1:M
%             XA_VAL=[XA_VAL;A(1,i+j)-A(1,i)];
%             YA_VAL=[YA_VAL;A(2,i+j)-A(2,i)];
%             ZA_VAL=[ZA_VAL;A(3,i+j)-A(3,i)];
%        end
         for i_m=1:M
              A_VAL(A_index,i_m)=A(i_m,i+j)-A(i_m,i);
         end
         A_index=A_index+1;
       end
    end
end

% x=A(1,:);
% a=repmat(x,[60,1]);
% result=a.*yBASE;
% DBASE=[base;2*base*t0];
% BASE=BASE(:,1:9980-1);
% BASE=BASE(:,1:10000-1);
% BASE=BASE(:,1:9988-1);
BASE_pinv=pinv(BASE);
% yBASE_pinv=pinv(yBASE);
for i_m=1:M
    diff(:,i_m)=BASE_pinv*A_VAL(:,i_m);
end
% diff_x=BASE_pinv*XA_VAL;
% diff_y=BASE_pinv*YA_VAL;
% diff_z=BASE_pinv*ZA_VAL;
%用得出的多项式对序列的预测
%  A_yu(1,1)=A(1,1);
%  A_yu(2,1)=A(2,1);
%  A_yu(3,1)=A(3,1);
 for i_m=1:M

    A_yu(i_m,1)=A(i_m,12000);
end
for i=1:2000
    base=baseEuqation(A_yu(:,i)')';
    yu_BASE=[base*t0,base*t0*t0];
%     diff_x*BASE
%     A_yu(1,i+1)=A_yu(1,i)+yu_BASE*diff_x;
%     A_yu(2,i+1)=A_yu(2,i)+yu_BASE*diff_y;
%     A_yu(3,i+1)=A_yu(3,i)+yu_BASE*diff_z;
    for i_m=1:M
%         diff(:,i_m)=BASE_pinv*A_VAL(:,i_m);  % 这行代码不需要!
        A_yu(i_m,i+1)=A_yu(i_m,i)+yu_BASE*diff(:,i_m);
    end
end
figure
aa_yu=[X_n(1,10000:12000),A_yu(1,1:2000)];
TTT=[1:4000];
plot(TTT,X_n(1,10000:14000-1),'LineWidth',3)
hold on
plot(TTT(1,2000:4000),aa_yu(1,2000:4000),'r --','LineWidth',3)
    set(gca,'FontSize',20);
xlabel('Step','fontsize',36);
ylabel('x','fontsize',36); 
set(gca,'linewidth',5)%加粗坐标轴
set(gca,'Xcolor',[0 0 0]);
set(gca,'Ycolor',[0 0 0]);


% save('ONLYZnoNoise_x_t0=0.01_L=10001.mat','A');
% save('ONLYZdiff_x.mat','diff_x');
% save('ONLYZdiff_y.mat','diff_y');
% save('ONLYZdiff_z.mat','diff_z');
%  A_yu(1,1)=A(1,1);
%  A_yu(2,1)=A(2,1);
%  A_yu(3,1)=A(3,1);
% for i=1:N
%     base=baseEuqation(A_yu(:,i)');
%     BASE=[base*t0;base*t0*t0];
% %     diff_x*BASE
%     A_yu(1,i+1)=A_yu(1,i)+diff_x*BASE;
%     A_yu(2,i+1)=A_yu(2,i)+diff_y*BASE;
%     A_yu(3,i+1)=A_yu(3,i)+diff_z*BASE;
% end
%在不知道动力学方程时，用多项式拟合
D_3=[];
D_3_m=[];
%   Dy_3=[];
%   Dz_3=[];
D_index=1;
for i=m+1:size(A,2)-m
    for j=-m:m
       D_base=baseEuqation(A(:,i)',m);
       D_BASE=[D_base,2*D_base*t0*j];
%        D_3=[D_3;D_BASE*diff_x,D_BASE*diff_y,D_BASE*diff_z];
       for i_m=1:M
           D_3_m(D_index,i_m)=D_BASE*diff(:,i_m);
       end
       D_index=D_index+1;
%        Dy_3=[Dy_3;diff_y*D_BASE];
%        Dz_3=[Dz_3;diff_z*D_BASE];
    end
end
% 
 Xbase=Z2baseEuqation(A',m);
% %  Ybase=ZbaseEuqation(A');
% %  Zbase=ZbaseEuqation(A');
% 
% 
  XBASE_pinv=pinv(Xbase);
%   YBASE_pinv=pinv(Ybase);
%   ZBASE_pinv=pinv(Zbase);
%   Dx_3=[];
%   Dy_3=[];
%   Dz_3=[];
%   m=1;
% A_f_x=load('diff_x');
% A_f_y=load('diff_y');
% A_f_z=load('diff_z');
% a_f_x=A_f_x.diff_x;
% a_f_y=A_f_y.diff_y;
% a_f_z=A_f_z.diff_z;
% % a_f=a_l.c_f;
% % N=L;
% syms t x y z
% % X=@(x,y,z,t)x+(a_f(1)*x+a_f(2)*y+a_f(3)*z+a_f(4)*x*x+a_f(5)*x*y+a_f(6)*x*z+a_f(7)*y*y+a_f(8)*y*z+a_f(9)*z*z+a_f(10)*x*x*x+a_f(11)*x*x*y+a_f(12)*x*x*z++a_f(13)*x*y*y+a_f(14)*x*z*z+a_f(15)*x*y*z+a_f(16)*y*y*y+a_f(17)*y*y*z+a_f(18)*z*z*z+a_f(19)*z*z*y)*t+(a_f(20)*x+a_f(21)*y+a_f(22)*z+a_f(23)*x*x+a_f(24)*x*y+a_f(25)*x*z+a_f(26)*y*y+a_f(27)*y*z+a_f(28)*z*z+a_f(29)*x*x*x+a_f(30)*x*x*y+a_f(31)*x*x*z++a_f(32)*x*y*y+a_f(33)*x*z*z+a_f(34)*x*y*z+a_f(35)*y*y*y+a_f(36)*y*y*z+a_f(37)*z*z*z+a_f(38)*z*z*y)*t*t;
% % Y=@(x,y,z,t)y+(a_f(39)*x+a_f(40)*y+a_f(41)*z+a_f(42)*x*x+a_f(43)*x*y+a_f(44)*x*z+a_f(45)*y*y+a_f(46)*y*z+a_f(47)*z*z+a_f(48)*x*x*x+a_f(49)*x*x*y+a_f(50)*x*x*z++a_f(51)*x*y*y+a_f(52)*x*z*z+a_f(53)*x*y*z+a_f(54)*y*y*y+a_f(55)*y*y*z+a_f(56)*z*z*z+a_f(57)*z*z*y)*t+(a_f(58)*x+a_f(59)*y+a_f(60)*z+a_f(61)*x*x+a_f(62)*x*y+a_f(63)*x*z+a_f(64)*y*y+a_f(65)*y*z+a_f(66)*z*z+a_f(67)*x*x*x+a_f(68)*x*x*y+a_f(69)*x*x*z++a_f(70)*x*y*y+a_f(71)*x*z*z+a_f(72)*x*y*z+a_f(73)*y*y*y+a_f(74)*y*y*z+a_f(75)*z*z*z+a_f(76)*z*z*y)*t*t;
% % Z=@(x,y,z,t)z+(a_f(77)*x+a_f(78)*y+a_f(79)*z+a_f(80)*x*x+a_f(81)*x*y+a_f(82)*x*z+a_f(83)*y*y+a_f(84)*y*z+a_f(85)*z*z+a_f(86)*x*x*x+a_f(87)*x*x*y+a_f(88)*x*x*z++a_f(89)*x*y*y+a_f(90)*x*z*z+a_f(91)*x*y*z+a_f(92)*y*y*y+a_f(93)*y*y*z+a_f(94)*z*z*z+a_f(95)*z*z*y)*t+(a_f(96)*x+a_f(97)*y+a_f(98)*z+a_f(99)*x*x+a_f(100)*x*y+a_f(101)*x*z+a_f(102)*y*y+a_f(103)*y*z+a_f(104)*z*z+a_f(105)*x*x*x+a_f(106)*x*x*y+a_f(107)*x*x*z++a_f(108)*x*y*y+a_f(109)*x*z*z+a_f(110)*x*y*z+a_f(111)*y*y*y+a_f(112)*y*y*z+a_f(113)*z*z*z+a_f(114)*z*z*y)*t*t;
% X=@(x,y,z,t)x+(a_f_x(1)+a_f_x(2)*x+a_f_x(3)*y+a_f_x(4)*z+a_f_x(5)*x*x+a_f_x(6)*x*y+a_f_x(7)*x*z+a_f_x(8)*y*y+a_f_x(9)*y*z+a_f_x(10)*z*z+a_f_x(11)*x*x*x+a_f_x(12)*x*x*y+a_f_x(13)*x*x*z++a_f_x(14)*x*y*y+a_f_x(15)*x*z*z+a_f_x(16)*x*y*z+a_f_x(17)*y*y*y+a_f_x(18)*y*y*z+a_f_x(19)*z*z*z+a_f_x(20)*z*z*y)*t+(a_f_x(21)+a_f_x(22)*x+a_f_x(23)*y+a_f_x(24)*z+a_f_x(25)*x*x+a_f_x(26)*x*y+a_f_x(27)*x*z+a_f_x(28)*y*y+a_f_x(29)*y*z+a_f_x(30)*z*z+a_f_x(31)*x*x*x+a_f_x(32)*x*x*y+a_f_x(33)*x*x*z+a_f_x(34)*x*y*y+a_f_x(35)*x*z*z+a_f_x(36)*x*y*z+a_f_x(37)*y*y*y+a_f_x(38)*y*y*z+a_f_x(39)*z*z*z+a_f_x(40)*z*z*y)*t*t;
% Y=@(x,y,z,t)y+(a_f_y(1)+a_f_y(2)*x+a_f_y(3)*y+a_f_y(4)*z+a_f_y(5)*x*x+a_f_y(6)*x*y+a_f_y(7)*x*z+a_f_y(8)*y*y+a_f_y(9)*y*z+a_f_y(10)*z*z+a_f_y(11)*x*x*x+a_f_y(12)*x*x*y+a_f_y(13)*x*x*z++a_f_y(14)*x*y*y+a_f_y(15)*x*z*z+a_f_y(16)*x*y*z+a_f_y(17)*y*y*y+a_f_y(18)*y*y*z+a_f_y(19)*z*z*z+a_f_y(20)*z*z*y)*t+(a_f_y(21)+a_f_y(22)*x+a_f_y(23)*y+a_f_y(24)*z+a_f_y(25)*x*x+a_f_y(26)*x*y+a_f_y(27)*x*z+a_f_y(28)*y*y+a_f_y(29)*y*z+a_f_y(30)*z*z+a_f_y(31)*x*x*x+a_f_y(32)*x*x*y+a_f_y(33)*x*x*z+a_f_y(34)*x*y*y+a_f_y(35)*x*z*z+a_f_y(36)*x*y*z+a_f_y(37)*y*y*y+a_f_y(38)*y*y*z+a_f_y(39)*z*z*z+a_f_y(40)*z*z*y)*t*t;
% Z=@(x,y,z,t)z+(a_f_z(1)+a_f_z(2)*x+a_f_z(3)*y+a_f_z(4)*z+a_f_z(5)*x*x+a_f_z(6)*x*y+a_f_z(7)*x*z+a_f_z(8)*y*y+a_f_z(9)*y*z+a_f_z(10)*z*z+a_f_z(11)*x*x*x+a_f_z(12)*x*x*y+a_f_z(13)*x*x*z++a_f_z(14)*x*y*y+a_f_z(15)*x*z*z+a_f_z(16)*x*y*z+a_f_z(17)*y*y*y+a_f_z(18)*y*y*z+a_f_z(19)*z*z*z+a_f_z(20)*z*z*y)*t+(a_f_z(21)+a_f_z(22)*x+a_f_z(23)*y+a_f_z(24)*z+a_f_z(25)*x*x+a_f_z(26)*x*y+a_f_z(27)*x*z+a_f_z(28)*y*y+a_f_z(29)*y*z+a_f_z(30)*z*z+a_f_z(31)*x*x*x+a_f_z(32)*x*x*y+a_f_z(33)*x*x*z+a_f_z(34)*x*y*y+a_f_z(35)*x*z*z+a_f_z(36)*x*y*z+a_f_z(37)*y*y*y+a_f_z(38)*y*y*z+a_f_z(39)*z*z*z+a_f_z(40)*z*z*y)*t*t;
% dx=diff(X(x,y,z,t),t);%得到X的导数
% dy=diff(Y(x,y,z,t),t);
% dz=diff(Z(x,y,z,t),t);
%   for i=m+1:N-m
%       for j=-m:m
%        Dx_3=[Dx_3;subs(dx,[x,y,z,t],[A(1,i),A(2,i),A(3,i),t0*j])];
%        Dy_3=[Dy_3;subs(dy,[x,y,z,t],[A(1,i),A(2,i),A(3,i),t0*j])];
%        Dz_3=[Dz_3;subs(dz,[x,y,z,t],[A(1,i),A(2,i),A(3,i),t0*j])];
% %   Dx_3=[Dx_3;Dx(1,m);Dx(1,m);Dx(1,m);];
% %   Dy_3=[Dy_3;Dy(1,m);Dy(1,m);Dy(1,m);];
% %   Dz_3=[Dz_3;Dz(1,m);Dz(1,m);Dz(1,m);];
% 
%       end
%   end
  D_diff=XBASE_pinv*D_3_m;
%   XYZdiff=XBASE_pinv*D_3;
%   Xdiff=XYZdiff(:,1)';
%   Ydiff=XYZdiff(:,2)';
%   Zdiff=XYZdiff(:,3)';

  %用构造的方程来算dx,dy,dz,并和D_3中的数据对比

xu_A=Xbase*D_diff(:,2);


%   Ydiff=D_3'*YBASE_pinv;
%   Zdiff=D_3'*ZBASE_pinv;
%   XYZbase=AbaseEuqationm_1(A);
%     X_VALUE=Xdiff* XYZbase;
%     Y_VALUE=Ydiff* XYZbase; 
%     Z_VALUE=Zdiff* XYZbase;
%     zbiao=1;
%     for i=2:3:size(D_3,1)
%         Dz(1,zbiao)=D_3(i,3);
%         zbiao=zbiao+1;
%     end
    
% 误差计算
% (sum(Dz.*Dz,2)-sum(Z_VALUE.*Z_VALUE,2))/size(Dz,2)
% m_err=Dz-Z_VALUE(2:size(A,2)-1);
% % % sum(m_err.*m_err,2)
% % % size(m_err)
% % sum(abs(m_err))
% sum(abs(m_err))/size(m_err,2)
% 
% err=sum(m_err.*m_err,2)/size(m_err,2);

%  
%解方程的固定点
% a_l=load('ONLYZc_f_m=1_t0=0.01_N=9988.mat','c_f');
% A_L=a_l.c_f;
% Xdiff=A_L(1,1:10);
% Ydiff=A_L(1,11:20);
% Zdiff=A_L(1,21:30);

% x00 = [0.1;0.1;0.1]; 
x00 =[];
exitflag=0;
% while exitflag==0
    x00 =[];
    for i=1:M
        
        x00=[x00;0];
%     x00=[x00;9];
    end% Make a starting guess at the solution
    options = optimoptions('fsolve','Display','iter'); % Option to display output
% [x_fix,fval] = fsolve(@(x)nolinerfun(x,Xdiff,Ydiff,Zdiff),x00,options); % Call solver
    [x_fix,fval,exitflag] = fsolve(@(x)nolinerfun(x,D_diff,M),x00,options); % Call solver
%     [x_fix,fval,exitflag] = fsolve(@(x)nolinerfun(x,D_diff,M),x00); % Call solver


% end
% A_M=[Xdiff(1:3);
%     Ydiff(1:3);
%     Zdiff(1:3)];
% trac=trace(A_M);
% D_et=det(A_M);
%在固定点出求雅可比矩阵
B_G=[];
jbbb=jaco(M,D_diff,x_fix);
jbbb=double(jbbb);
[X_eig,B_eig]=eig(jbbb);
for i=1:M
%      if (isreal(B_eig(i,i)))
         B_G=[B_G;B_eig(i,i)];
%      end
end
% B_G=sort(B_G,'descend');
% trac=sum(B_G(1:3,1));
% D_et=B_G(1,1)*B_G(2,1)*B_G(3,1);
% % trac=double(trace(jbbb));
% % D_et=double(det(jbbb));
% 
% % B_M=[Xdiff(1,2)+2*Xdiff(1,5)*x_fix(1,1)+Xdiff(1,6)*x_fix(2,1)+Xdiff(1,7)*x_fix(3,1),Xdiff(1,3)+Xdiff(1,6)*x_fix(1,1)+2*Xdiff(1,8)*x_fix(2,1)+Xdiff(1,9)*x_fix(3,1),Xdiff(1,4)+Xdiff(1,7)*x_fix(1,1)+Xdiff(1,9)*x_fix(2,1)+2*Xdiff(1,10)*x_fix(3,1);
% % Ydiff(1,2)+2*Ydiff(1,5)*x_fix(1,1)+Ydiff(1,6)*x_fix(2,1)+Ydiff(1,7)*x_fix(3,1),Ydiff(1,3)+Ydiff(1,6)*x_fix(1,1)+2*Ydiff(1,8)*x_fix(2,1)+Ydiff(1,9)*x_fix(3,1),Ydiff(1,4)+Ydiff(1,7)*x_fix(1,1)+Ydiff(1,9)*x_fix(2,1)+2*Ydiff(1,10)*x_fix(3,1);
% % Zdiff(1,2)+2*Zdiff(1,5)*x_fix(1,1)+Zdiff(1,6)*x_fix(2,1)+Zdiff(1,7)*x_fix(3,1),Zdiff(1,3)+Zdiff(1,6)*x_fix(1,1)+2*Zdiff(1,8)*x_fix(2,1)+Zdiff(1,9)*x_fix(3,1),Zdiff(1,4)+Zdiff(1,7)*x_fix(1,1)+Zdiff(1,9)*x_fix(2,1)+2*Zdiff(1,10)*x_fix(3,1)];
% % trac=trace(B_M);
% % D_et=det(B_M);
% %根据相似矩阵的性质，求方程的参数
% syms x_x y_y
% [r1,r2]=solve(27*x_x*y_y==D_et,-x_x-1-y_y==trac,x_x,y_y);%洛伦兹方程
% % [r1,r2]=solve(xy==D,x-y==t1,x,y);
% reult=[r1,r2];
% vpa(reult)
% 
syms a_a c_c
%若斯勒方程
[r1,r2]=solve(a_a-a_a*x_fix(1,1)-c_c==-5.69015524977061+0.0971585415918918 + 0.995055709481755i+0.0971585415918918 - 0.995055709481755i,-a_a*x_fix(1,1)-c_c-x_fix(1,1)*a_a==-5.69015524977061*(0.0971585415918918 + 0.995055709481755i)*(0.0971585415918918 - 0.995055709481755i),a_a,c_c);%若斯勒方程
% [r1,r2]=solve(xy==D,x-y==t1,x,y);
reult=[r1,r2];
vpa(reult)
b_b=-x_fix(1,1)*reult(1)*x_fix(1,1)-x_fix(1,1)*reult(2);
vpa(b_b)
% 
% 
% 
% bx=a*A(1,:);
% DDZ=Dx+bx;
% diff_yy=DDZ*yBASE_pinv;
% % X=A(1,1:end-1)+diff_x*BASE;%得到拟合后的x分量的数据
% % Y=A(2,1:end-1)+diff_y*BASE;
% % Z=A(3,1:end-1)+diff_z*BASE;
% % result=[result;-A(3,:)];
% % zBASE_pinv=pinv(result);
% % diff_yy=Dz*zBASE_pinv;
% 
% yBASE=[base;base*t0;base*t0*t0];
% YY=diff_yy*yBASE;
% 
% Dy=diff_yy(21:60)*DBASE;
% XYY=A(2,:)-A(1,:);
% XYY_pinv=pinv(XYY);
% A_value=Dx*XYY_pinv;
% % BKK=YY-A(1,:);
% % BKK_PI=pinv(BKK);
% zBASE=[-base;-base*t0;-base*t0*t0];
% x=A(1,:);
% a=repmat(x,[60,1]);
% result=a.*zBASE;
% result=[A(1,:);result];
% DDY=YY+Dy;
% ZBASE_pinv=pinv(result);
% zz_diff=DDY*ZBASE_pinv;
% 
% plot3(X,YY(1,1:end-1),Z);
% plot3(X,Y,Z);