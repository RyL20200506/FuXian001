%%
clear all;clc;
% Lorenz system
x_0=-8;y_0=7;z_0=27;
x_c=10;y_c=8/3;z_c=28;
dt_1=0.01;n=5000;N=1;m=1;tau=0;nn=8000;
[x,y,z]=Runge_Kutta(dt_1,nn,x_0,y_0,z_0,x_c,y_c,z_c);
%取时序中N:n这一段做拟合
L1(1,:)=x(N:n);
L1(2,:)=y(N:n);
L1(3,:)=z(N:n);
L(1,:)=L1(1,1:end-1);L(2,:)=L1(2,1:end-1);L(3,:)=L1(3,1:end-1);
L41(1,:)=L1(1,2:end);L41(2,:)=L1(2,2:end);L41(3,:)=L1(3,2:end);
LL(1,:)=x(n:n+3000);LL(2,:)=y(n:n+3000);LL(3,:)=z(n:n+3000);

var=baseEuqation1(L');% g_i,f_i
var_1=[var(:,1:end).*(dt_1);var(:,1:end).*(dt_1).^2];
%
B_1=[L41(1,1:end)-L(1,1:end);L41(2,1:end)-L(2,1:end);L41(3,1:end)-L(3,1:end);]*pinv(var_1(:,1:end)); %用前n个数据估计系数B_1

subplot(2,1,1);
x_fitting=B_1(1,:)*[var(:,1:end).*(dt_1);var(:,1:end).*(dt_1).^2]+L(1,1:end);
z_fitting=B_1(3,:)*[var(:,1:end).*(dt_1);var(:,1:end).*(dt_1).^2]+L(3,1:end);
plot(x_fitting,'b');hold on;
plot(L41(1,1:end),'r');
legend('fitting data','real data');
xlabel('t');ylabel('x');
subplot(2,1,2);
plot(abs(x_fitting-L41(1,1:end)));hold on;
xlabel('t');ylabel('absolute error');

%生成新时序
figure;
x_appro=zeros(3,length(LL));
x_appro(1,1)=LL(1,1);%用第n+1个生成新的时序
x_appro(2,1)=LL(2,1);
x_appro(3,1)=LL(3,1);
A=zeros(length(x_appro),3);
for i=1:length(x_appro)-1
    
        A(i,1)=x_appro(1,i);
        A(i,2)=x_appro(2,i);
        A(i,3)=x_appro(3,i);
        var_0=baseEuqation1(A(i,:));
        var_00=[var_0*dt_1;var_0*dt_1.^2];
        x_appro(:,i+1)=B_1*var_00+(x_appro(:,i));

end
% 新的时序与原本时序LL进行画图对比
t_1=[1:(length(x_appro))];
plot(x_appro(1,1:end),'b');hold on;
t_2=[1:(length(L)+length(LL))];
plot(LL(1,1:end),'r--','Linewidth',1.1)
legend('Predict','Lorenz');
xlabel('step');ylabel('x');