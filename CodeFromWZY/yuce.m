close all;clear;clc
for h=1:10
a=10;b=8/3;r=28;
x0=[-8;7;27];
xn=[x0];
D=0;t_u=0.01;N=10000;
b1=sqrt(2*D*t_u)*t_u;b2=0.5*sqrt(2*D*t_u)*t_u;b3=(1/2*sqrt(3))*sqrt(2*D*t_u)*t_u;
w1=randn(3,N);w2=randn(3,N);
%--------------------------根据公式A3求解运动轨迹-------------------------------------
for i=1:N
    t1=b1*w1(:,i);
    t2=b2*w1(:,i)+b3*w2(:,i);
    dx=a*(xn(2,i)-xn(1,i));
    dy=r*xn(1,i)-xn(2,i)-xn(1,i)*xn(3,i);
    dz=xn(1,i)*xn(2,i)-b*xn(3,i);
    df1=[a*(dy-dx);r*dx-dy-dx*xn(3,i)-xn(1,i)*dz;dx*xn(2,i)+xn(1,i)*dy-b*dz];
    df2=[a;r-xn(3,i)-1-xn(1,i);xn(2,i)+xn(1,i)-b];
    xt=xn(:,i)+[dx;dy;dz]*t_u+0.5*t_u*t_u*df1+df2.*t2+t1;
    xn=[xn,xt];
end
A=xn(:,1:end-1);
%subplot(1,2,1)
%plot(A(1,:),A(3,:),'linewidth',1);
%xlabel('x');
%ylabel('z');
%set(gca,'FontSize',30) ;
%subplot(1,2,2)
%plot(A(1,:),'linewidth',1);
%xlabel('step');
%ylabel('x');
%set(gca,'FontSize',30) ;
B=A(:,1:5000);
%------------------------------使用公式3对g，h进行拟合----------------------------------
BASE=[];
base=baseEuqation(B');%g、h系数后的多项式
BASE=[base*t_u;base*t_u*t_u];%m=1,从第一个元素算起,初始点数为1：L，对2：L+1拟合的结果
DBASE=[base;2*t_u*base];%对t求导
BASE=BASE(:,1:end-1);
BASE_pinv=pinv(BASE);
%---------------------------------多项式前的系数---------------------------------------
cof_x=(B(1,2:end)-B(1,1:end-1))*BASE_pinv;
cof_y=(B(2,2:end)-B(2,1:end-1))*BASE_pinv;
cof_z=(B(3,2:end)-B(3,1:end-1))*BASE_pinv;
%---------------------------------下一个点的值-----------------------------------------
X=B(1,1:end-1)+cof_x*BASE;
Y=B(2,1:end-1)+cof_y*BASE;
Z=B(3,1:end-1)+cof_z*BASE;
%---------------------------------对x,y,z进行求导--------------------------------------
fx=cof_x*DBASE;
fy=cof_y*DBASE;
fz=cof_z*DBASE;
%----------------------------洛伦兹方程求a,r,b相乘的两个矩阵----------------------------
for i=1:5000
    base_a(i,:)=B(2,i)-B(1,i);
end
for i=1:5000
    base_r(i,:)=B(1,i);
end
for i=1:5000
    base_b(i,:)=B(3,i);
end
for i=1:5000
    pol_a(i,:)=fx(1,i);
end
for i=1:5000
    pol_r(i,:)=fy(1,i)+B(2,i)+B(1,i)*B(3,i);
end
for i=1:5000
    pol_b(i,:)=B(1,i)*B(2,i)-fz(1,i);
end
pinv_a=pinv(base_a);
pinv_r=pinv(base_r);
pinv_b=pinv(base_b);
%---------------------------------------求值--------------------------------------------
a_pre(h,1)=pinv_a*pol_a;
r_pre(h,1)=pinv_r*pol_r;
b_pre(h,1)=pinv_b*pol_b;
for i=1:size(B,2)%i=1:矩阵A的列数
    base_a(i,:)=(B(2,i)-B(1,i))*t_u;
end
 base_ap=pinv(base_a);%x'
 BASE_11=[base*t_u;base*t_u*t_u];%系数后多项式
    for i=1:size(B,2)
       DERA(i,:)=0.5*t_u*t_u*(-a_pre(h,1)*a_pre(h,1)*(B(2,i)-B(1,i))+a_pre(h,1)*(r_pre(h,1)*B(1,i)-B(2,i)-B(1,i)*B(3,i)));
    %0.5*t0*t0*(a*(y'-x'))=0.5*t0*t0*x''
    end
     a_pre(h,1)=base_ap*(BASE_11'*cof_x'-DERA);%二阶泰勒展开%a=[(y-x)*t0]的逆*[x(t+t0)-x(t)-0.5*t0*t0*x'']
    for i=1:size(B,2)
        base_y(i,:)=B(1,i)*t_u;
    end
    base_yp=pinv(base_y);
    for i=1:size(B,2)
       DERA_y(i,:)=0.5*t_u*t_u*( a_pre(h,1)*(B(2,i)-B(1,i))*(r_pre(h,1)-B(3,i))-( r_pre(h,1)*B(1,i)-B(2,i)-B(1,i)*B(3,i))-B(1,i)*(B(2,i)*B(1,i)-b_pre(h,1)*B(3,i)));
    end
    
    r_pre(h,1)=base_yp*(BASE_11'*cof_y'-DERA_y+t_u*(B(2,:)+B(1,:).*B(3,:))');
    for i=1:size(B,2)
        base_z(i,:)=B(3,i)*t_u;
    end
    base_zp=pinv(base_z);
    for i=1:size(B,2)
       DERA_z(i,:)=0.5*t_u*t_u*(a_pre(h,1)*(B(2,i)-B(1,i))*B(2,i)+(r_pre(h,1)*B(1,i)-B(2,i)-B(1,i)*B(3,i))*A(1,i)-b_pre(h,1)*(B(2,i)*B(1,i)-b_pre(h,1)*B(3,i)));
    end
    b_pre(h,1)=-base_zp*(BASE_11'*cof_z'-DERA_z-t_u*(B(1,:).*B(2,:))');
end
a_value=mean(a_pre)
r_value=mean(r_pre)
b_value=mean(b_pre)
%----------------------------------预测轨迹---------------------------------------------
for i=1:N
    t1=b1*w1(:,i);
    t2=b2*w1(:,i)+b3*w2(:,i);
    dx=a_value*(xn(2,i)-xn(1,i));
    dy=r_value*xn(1,i)-xn(2,i)-xn(1,i)*xn(3,i);
    dz=xn(1,i)*xn(2,i)-b_value*xn(3,i);
    df1=[a_value*(dy-dx); r_value*dx-dy-dx*xn(3,i)-xn(1,i)*dz; dx*xn(2,i)+xn(1,i)*dy-b_value*dz];
    df2=[a_value; r_value-xn(3,i)-1-xn(1,i); xn(2,i)+xn(1,i)-b_value];
    xt=xn(:,i)+[dx;dy;dz]*t_u+0.5*t_u*t_u*df1+df2.*t2+t1;
    x_n=[xn,xt];
end
AA=x_n(:,1:end-1);
plot(AA(1,:),'b --','linewidth',2);
hold on;
plot(A(1,:),'r -','linewidth',1);
hold off;
xlabel('step');
ylabel('x');
set(gca,'FontSize',30) ;

%----------------------------------调试轨迹---------------------------------------------
for i=1:N
    xt=xn(:,i);
    x_n=[xn,xt];
end
AA=x_n(:,1:end-1);
plot(AA(1,:),'b --','linewidth',2);
hold on;
plot(A(1,:),'r -','linewidth',1);
hold off;
xlabel('step');
ylabel('x');
set(gca,'FontSize',30) ;