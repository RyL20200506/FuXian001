%% Taylor cacluator nosie Lorenz System
close all;clear;clc;

delta=10;b=8/3;r=28;
N=15000;
x0=[-8;7;27];
X_n=[x0];
t0=0.01;%单位计算时间
m=1;
D=0;%噪声方差
% D=0;


% 正常得到X_n是真实的时间序列
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

B=X_n(:,1:12000);

k_length=length(B);
%%%%重构维度的构建

M=3;
A=B;
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
         for i_m=1:M
              A_VAL(A_index,i_m)=A(i_m,i+j)-A(i_m,i);
         end
         A_index=A_index+1;
       end
    end
end

BASE_pinv=pinv(BASE);
for i_m=1:M
    diff(:,i_m)=BASE_pinv*A_VAL(:,i_m);  % 疑: 这个diff是不是使用了全部的12000条数据才拟合出来的? 那是不是使用了未来函数
end

 for i_m=1:M
    A_yu(i_m,1)=A(i_m,5000);  % 先给前5000列, 表示前5000个时刻点.
 end

for i=1:size(A,2)
    base=baseEuqation(A_yu(:,i)')';
    yu_BASE=[base*t0,base*t0*t0];
    for i_m=1:M
        diff(:,i_m)=BASE_pinv*A_VAL(:,i_m);  % 怀疑是这里!!!很严重的使用了未来函数
        A_yu(i_m,i+1)=A_yu(i_m,i)+yu_BASE*diff(:,i_m);
    end
end

aa_yu=[X_n(1,1:5000),A_yu(1,1:10000)];
TTT=[1:10000];
plot(TTT,X_n(1,1:10000),'LineWidth',3)
hold on
plot(TTT(1,5000:10000),aa_yu(1,5000:10000),'r --','LineWidth',3)
    set(gca,'FontSize',20);
xlabel('Step','fontsize',36);
ylabel('x','fontsize',36); 
set(gca,'linewidth',5)%加粗坐标轴
set(gca,'Xcolor',[0 0 0]);
set(gca,'Ycolor',[0 0 0]);


D_3=[];
D_3_m=[];
D_index=1;
for i=m+1:size(A,2)-m
    for j=-m:m
       D_base=baseEuqation(A(:,i)',m);
       D_BASE=[D_base,2*D_base*t0*j];
       for i_m=1:M
           D_3_m(D_index,i_m)=D_BASE*diff(:,i_m);
       end
       D_index=D_index+1;
    end
end

% 
Xbase=Z2baseEuqation(A',m);
XBASE_pinv=pinv(Xbase);
D_diff=XBASE_pinv*D_3_m;
xu_A=Xbase*D_diff(:,2);

x00 =[];
exitflag=0;
    x00 =[];
    for i=1:M
        x00=[x00;0];
    end% Make a starting guess at the solution
    options = optimoptions('fsolve','Display','iter'); % Option to display output
    [x_fix,fval,exitflag] = fsolve(@(x)nolinerfun(x,D_diff,M),x00,options); % Call solver

    
%在固定点出求雅可比矩阵
B_G=[];
jbbb=jaco(M,D_diff,x_fix);
jbbb=double(jbbb);
[X_eig,B_eig]=eig(jbbb);
for i=1:M
         B_G=[B_G;B_eig(i,i)];
end


syms a_a c_c
%若斯勒方程
[r1,r2]=solve(a_a-a_a*x_fix(1,1)-c_c==-5.69015524977061+0.0971585415918918 + 0.995055709481755i+0.0971585415918918 - 0.995055709481755i,-a_a*x_fix(1,1)-c_c-x_fix(1,1)*a_a==-5.69015524977061*(0.0971585415918918 + 0.995055709481755i)*(0.0971585415918918 - 0.995055709481755i),a_a,c_c);%若斯勒方程
% [r1,r2]=solve(xy==D,x-y==t1,x,y);
reult=[r1,r2];
vpa(reult)
b_b=-x_fix(1,1)*reult(1)*x_fix(1,1)-x_fix(1,1)*reult(2);
vpa(b_b)

