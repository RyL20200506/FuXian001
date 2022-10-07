close all;clear all;clc;
clear VAR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%榫欐牸搴撳
% t0=0.01;
% [t1,h1]=runge_kutta(@lor_fun,[12 4 8],t0,0,5000);
% X_n=h1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%浜岄樁闅忔満榫欐牸搴撳
close all;clear;clc;
% for r_value=26:2:64
for r_value=28  % 生成少一点的r
    delta=10;b=8/3;r=r_value;
    N=40000;
    % x0=[-1;3;4];
    x0=[12;4;8];
    X_n=[x0];
    t0=0.01;  % 鍗曚綅璁＄畻鏃堕棿
    m=1;
    % D=0.5;%鍣０鏂瑰樊
    D=0;
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

    NK=N;
    x_ce=X_n(1,500:NK)';
    x_1=X_n(1,499:NK)';
    YY=X_n(2,500:NK);
    ZZ=X_n(3,500:NK);
    % plot(x_ce,'LineWidth',1.3);
    % xlabel('step');
    % ylabel('x');
    % % title('kuramoto oscillator network',"FontName","Times New Roman","FontSize",35,'FontWeight','bold');
    % set(gca,'fontsize',20,'fontname','Times New Roman','FontWeight','bold');
    % set(gca,'FontSize',20); 
    % set(gca,'linewidth',3)%加粗坐标轴
    %  set(gca,'XLim',[10000 15000]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%准确计算
    % x_dot_q=15*(X_n(2,:)-X_n(1,:));
    % % % y_dot_q=35*X_n(1,:)-X_n(2,:)-X_n(1,:).*X_n(3,:);
    % % % x_dot2_q=10*(y_dot_q-x_dot_q);
    % x_1dot=x_dot_q(1,500:NK);
    % x_2dot=x_dot2_q(1,500:NK);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%王方法
    % NKK=490000;
    % % NKK=4;
    % [x_2,x_2dot,yu]=solve_dot(x_1,t0,NKK);
    % % % % x_2dot=solve_dot(x_1,t0)';
    % % x_dot2=x_2dot';
    % xx_1=yu(1,499:NKK)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
    % x_ce=[];
    % x_ce=yu(1,500:NKK)';
    % % for i=1:NKK-1000
    % %     x_dot(i,:)=(xx_1(i+2,:)-xx_1(i,:))/(2*t0);
    % % end                            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % x_2dot=[];
    % [x_2,x_2dot,yu]=solve_dot(x_1,t0,1);
    % x_dot2=x_2dot';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 
    for i=1:NK-1000
        x_dot(i,:)=(x_1(i+2,:)-x_1(i,:))/(2*t0);
    end
    % x_1=[];
    % x_1=x_2';
    % x_dot=x_1dot';
    % x_2dot=x_2dot';
    %%%%%%
    % xita(1:NK-1000,1)=15;
    xita(1,1)=1;
    % rr(1:NK-1000,1)=35;
    rr(1,1)=1;
    % bb(1:NK-1000,1)=3;
    bb(1,1)=1;
    a1(1,1)=0;
    r1(1,1)=0;
    b1(1,1)=0;
    beta1=15;
    beta2=15;
    beta3=15;
    % y=zeros(1,15000);
    % z=zeros(1,15000);
    strr=num2str(r_value);
    wenjianname=strcat('x_noise0_r',strr,'.mat');
    save(wenjianname,'x_ce');
    % save('x2_noise5.mat','xx_1');
    wenjianname1=strcat('x_dot_noise0_r',strr,'.mat');
    save(wenjianname1,'x_dot');
    % save('x_dot2_noise5.mat','x_dot2');
    wenjiannamey=strcat('Y_noise0_r',strr,'.mat');
    wenjiannameZ=strcat('Z_noise0_r',strr,'.mat');
    save(wenjiannamey,'YY');
    save(wenjiannameZ,'ZZ');
end
x_ceq=load('x.mat');
x_dotq=load('x_dot.mat');
x_ce=x_ceq.x_ce;
x_dot=x_dotq.x_dot;
z(1,1)=1;
y(1,1)=1;
gggg=[];
ggggg=[];
gggggg=[];
g7=[];
g8=[];
k=1;
g_xita=[];
theta=1;
a=delta;
global VAL;

 [t1,H1]=mod_runge_kutta(@new_lor_lorenz,[0 0 1 1 1 0 1 0],t0,0,120,x_ce,x_dot,a);
% [t1,H1]=mod_runge_kutta(@mod_lor_fun,[0 0 1 1 1 1 1 0 0 0 0.1],t0,0,120,x_ce,x_dot);
for i=1:6
    plot(H1(i,1:60));
    hold on
end
for i=1:10000
    g(i,1)=x_dot(i,1)-H1(5,i)*(H1(3,i)-x_ce(i,1));
    y_err(i,1)=abs(H1(3,i)-YY(1,i));
end
figure
 plot(H1(3,5000:10000));
 hold on
 plot(YY(5000:10000));
 figure
 plot(H1(4,1:10000));
 hold on
 plot(ZZ(1:10000));
 
 figure
 plot(H1(4,:),'LineWidth',2);
 hold on
 plot(H1(3,:),'LineWidth',2);
 hold on
 plot(H1(5,:),'LineWidth',2);
  set(gca,'FontSize',20);
xlabel('Step','fontsize',36);
ylabel('value','fontsize',36); 
set(gca,'linewidth',5)%鍔犵矖鍧愭爣杞?
set(gca,'Xcolor',[0 0 0]);
set(gca,'Ycolor',[0 0 0]);
 figure
 subplot(1,2,1);
 plot(H1(1,3000:5000),'LineWidth',3);
   set(gca,'FontSize',20);
xlabel('Step','fontsize',36);
ylabel('y_{value}','fontsize',36); 
set(gca,'linewidth',5)%鍔犵矖鍧愭爣杞?
set(gca,'Xcolor',[0 0 0]);
set(gca,'Ycolor',[0 0 0]);
 hold on
 plot(YY(1,3000:5000),'LineWidth',1);
 subplot(1,2,2);
%  plot(H1(2,3000:5000),'LineWidth',3);
%  hold on
%  plot(ZZ(1,3000:5000),'LineWidth',1);
    set(gca,'FontSize',20);
xlabel('Step','fontsize',36);
ylabel('z_{value}','fontsize',36); 
set(gca,'linewidth',5)%鍔犵矖鍧愭爣杞?
set(gca,'Xcolor',[0 0 0]);
set(gca,'Ycolor',[0 0 0]);
 
 
 plot(g(1:5000,1),'b','LineWidth',2)
set(gca,'FontSize',20);
xlabel('Step','fontsize',36);
ylabel('absolute error','fontsize',36); 
set(gca,'linewidth',5)%鍔犵矖鍧愭爣杞?
set(gca,'Xcolor',[0 0 0]);
set(gca,'Ycolor',[0 0 0]);
title('Error');
ylim([-200 200])

mean(H1(4,6000:end))
mean(H1(3,6000:end))
mean(H1(5,6000:end))



