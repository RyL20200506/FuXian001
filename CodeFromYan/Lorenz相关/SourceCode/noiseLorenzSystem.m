% Taylor cacluator nosie Lorenz System
close all;clear;clc;
% for D=1:3
for zkkk=1:10
delta=10;b=8/3;r=28;
N=10000;
% x0=[-1;3;4];
x0=[-8;7;27];
X_n=[x0];
t0=0.01;%单位计算时间

D=50;%噪声方差
m=1;
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
A=X_n(:,1:end-1);
subplot(1,2,1)
plot(A(1,:),A(3,:),'linewidth',1);
xlabel('x');
ylabel('z');
set(gca,'FontSize',30) ;
subplot(1,2,2)
plot(A(1,:),'linewidth',1);
xlabel('step');
ylabel('x');
set(gca,'FontSize',30) ;
BASE=[];
base=baseEuqation(A');
BASE=[base*t0;base*t0*t0];%m=1,从第一个元素算起,初始点数为1：L，对2：L+1拟合的结果
DBASE=[base;2*t0*base];
BASE=BASE(:,1:N-1);
BASE_pinv=pinv(BASE);
diff_x=(A(1,2:end)-A(1,1:end-1))*BASE_pinv;
diff_y=(A(2,2:end)-A(2,1:end-1))*BASE_pinv;
diff_z=(A(3,2:end)-A(3,1:end-1))*BASE_pinv;
X=A(1,1:end-1)+diff_x*BASE;%得到拟合后的x分量的数据
% X=diff_x*BASE;
DX=diff_x*DBASE;
DY=diff_y*DBASE;
DZ=diff_z*DBASE;
Y=A(2,1:end-1)+diff_y*BASE;
% Y=diff_y*BASE;
Z=A(3,1:end-1)+diff_z*BASE;
% Z=diff_z*BASE;
% plot(X_n(1,:),X_n(2,:))
    for i=1:N
        base_a(i,:)=(A(2,i)-A(1,i));
    end
    for i=1:N
        base_b(i,:)=A(3,i);
    end
    for i=1:N
        base_r(i,:)=A(1,i);
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
    
    base_ap=pinv(base_a);
    base_bp=pinv(base_b);
    base_rp=pinv(base_r);
 for i=1:N
        cha_x(i,:)=DX(1,i);
 end
 for i=1:N
        cha_z(i,:)=A(1,i)*A(2,i)-DZ(1,i);
 end
  for i=1:N
        cha_y(i,:)=DY(1,i)+A(2,i)+A(1,i)*A(3,i);
 end
 a_value(zkkk,1)=base_ap*cha_x;
 b_value(zkkk,1)=base_bp*cha_z;
 r_value(zkkk,1)=base_rp*cha_y;
 %用论文11式精确方程
 for i=1:size(A,2)
    base_a(i,:)=(A(2,i)-A(1,i))*t0;
 end
 base_ap=pinv(base_a);
 BASE_11=[base*t0;base*t0*t0];
    for i=1:size(A,2)
       DERA(i,:)=0.5*t0*t0*(-a_value(zkkk,1)*a_value(zkkk,1)*(A(2,i)-A(1,i))+a_value(zkkk,1)*(r_value(zkkk,1)*A(1,i)-A(2,i)-A(1,i)*A(3,i)));
    end
     a_value(zkkk,1)=base_ap*(BASE_11'*diff_x'-DERA);
    for i=1:size(A,2)
        base_y(i,:)=A(1,i)*t0;
%         base_y(i,:)=[A(1,i)*t0,A(2,i)*t0,A(1,i)*A(3,i)*t0];
    end
    base_yp=pinv(base_y);
    for i=1:size(A,2)
       DERA_y(i,:)=0.5*t0*t0*( a_value(zkkk,1)*(A(2,i)-A(1,i))*(r_value(zkkk,1)-A(3,i))-( r_value(zkkk,1)*A(1,i)-A(2,i)-A(1,i)*A(3,i))-A(1,i)*(A(2,i)*A(1,i)-b_value(zkkk,1)*A(3,i)));
    end
    
    r_value(zkkk,1)=base_yp*(BASE_11'*diff_y'-DERA_y+t0*(A(2,:)+A(1,:).*A(3,:))');
    for i=1:size(A,2)
        base_z(i,:)=A(3,i)*t0;
    end
    base_zp=pinv(base_z);
    for i=1:size(A,2)
       DERA_z(i,:)=0.5*t0*t0*(a_value(zkkk,1)*(A(2,i)-A(1,i))*A(2,i)+(r_value(zkkk,1)*A(1,i)-A(2,i)-A(1,i)*A(3,i))*A(1,i)-b_value(zkkk,1)*(A(2,i)*A(1,i)-b_value(zkkk,1)*A(3,i)));
    end
    b_value(zkkk,1)=-base_zp*(BASE_11'*diff_z'-DERA_z-t0*(A(1,:).*A(2,:))');
 
 err_a(zkkk,1)=abs(a_value(zkkk,1)-delta)/delta;
 err_b(zkkk,1)=abs(b_value(zkkk,1)-b)/b;
 err_r(zkkk,1)=abs(r_value(zkkk,1)-r)/r;
 noise_DD=0;
  for i=1:N-1
       noise_D=(X(1,i)-A(1,i+1))^2;
       noise_DD= noise_DD+ noise_D;
  end
noise_DD= noise_DD/(N-1);
D_D(zkkk,1)=noise_DD/(2*t0);
err_D(zkkk,1)=abs(D_D(zkkk,1)-D)/D;
% delta=a_value;b=b_value;r=-r_value;
end
% a_err=(mean(a_value)-delta)/delta;
% b_err=(mean(b_value)-b)/b;
% r_err=(mean(r_value)-r)/r;
% D_err=(mean(D_D)-D)/D;
a_err=mean(err_a)
b_err=mean(err_b)
r_err=mean(err_r)
D_err=mean(err_D)

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