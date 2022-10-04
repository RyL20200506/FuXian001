% 学习第一部分, 是Lorenz的画图

% Taylor cacluator nosie Lorenz System
close all;clear;clc;  % R: 正常操作
% for D=1:3  
delta=10;b=8/3;r=28;  % Lorenz正常参数
N=10000;  
% x0=[-1;3;4];
x0 = [-8;7;27];  % 初始值
X_n = [x0];  % R: 容器
t0 = 0.01;  % 单位计算时间

D = 50;  % 噪声方差
m = 1;  % R: 哦这个好像是时间窗口之类的
b1=sqrt(2*D*t0);b2=1/2*b1*t0;b3=1/sqrt(3)*b2;  % ??? 这是啥 
w1=randn(N,3);w2=randn(N,3);  % R: 生成标准白噪声

% R: 下面应该是手写了数值积分.
for i=1:N
    % 
    s1=b1*w1(i,:);  % 
    s2=b2*w1(i,:)+b3*w2(i,:);  % 
    df2=[0;(r-X_n(3,i)-1-X_n(1,i));(X_n(2,i)+X_n(1,i)-b)];  %
    
    % R: 下面是Lorenz方程, 应该是原方程, 而不是拟合方程.
    fx=delta*(X_n(2,i)-X_n(1,i));  %
    fy=r*X_n(1,i)-X_n(2,i)-X_n(1,i)*X_n(3,i);  %
    fz=X_n(1,i)*X_n(2,i)-b*X_n(3,i);
    
    % R: 总之应该是往后推一段X_n
    df1=[-delta*fx+delta*fy; (r-X_n(3,i))*fx-fy-X_n(1,i)*fz; X_n(2,i)*fx+X_n(1,i)*fy-b*fz];  % 
    X_t=X_n(:,i)+[fx;fy;fz]*t0+1/2*t0*t0*df1+df2.*s2'+s1';  % R: 用最后一列数据去更新状态
    X_n=[X_n,X_t];  % R: 合并新的列
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
