% ѧϰ��һ����, ��Lorenz�Ļ�ͼ

% Taylor cacluator nosie Lorenz System
close all;clear;clc;  % R: ��������
% for D=1:3  
delta=10;b=8/3;r=28;  % Lorenz��������
N=10000;  
% x0=[-1;3;4];
x0 = [-8;7;27];  % ��ʼֵ
X_n = [x0];  % R: ����
t0 = 0.01;  % ��λ����ʱ��

D = 50;  % ��������
m = 1;  % R: Ŷ���������ʱ�䴰��֮���
b1=sqrt(2*D*t0);b2=1/2*b1*t0;b3=1/sqrt(3)*b2;  % ??? ����ɶ 
w1=randn(N,3);w2=randn(N,3);  % R: ���ɱ�׼������

% R: ����Ӧ������д����ֵ����.
for i=1:N
    % 
    s1=b1*w1(i,:);  % 
    s2=b2*w1(i,:)+b3*w2(i,:);  % 
    df2=[0;(r-X_n(3,i)-1-X_n(1,i));(X_n(2,i)+X_n(1,i)-b)];  %
    
    % R: ������Lorenz����, Ӧ����ԭ����, ��������Ϸ���.
    fx=delta*(X_n(2,i)-X_n(1,i));  %
    fy=r*X_n(1,i)-X_n(2,i)-X_n(1,i)*X_n(3,i);  %
    fz=X_n(1,i)*X_n(2,i)-b*X_n(3,i);
    
    % R: ��֮Ӧ����������һ��X_n
    df1=[-delta*fx+delta*fy; (r-X_n(3,i))*fx-fy-X_n(1,i)*fz; X_n(2,i)*fx+X_n(1,i)*fy-b*fz];  % 
    X_t=X_n(:,i)+[fx;fy;fz]*t0+1/2*t0*t0*df1+df2.*s2'+s1';  % R: �����һ������ȥ����״̬
    X_n=[X_n,X_t];  % R: �ϲ��µ���
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
