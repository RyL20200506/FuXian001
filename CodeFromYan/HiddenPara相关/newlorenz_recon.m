
close all;clear all;clc;
a_mean=[];
b_mean=[];
r_mean=[];
a_std=[];
b_std=[];
r_std=[];
index_value=1;
alfa=3;
learn_rate=0.005;
% alfa_value=1;
alfa=3;

t0=0.01;
x_ceq=load('x_noise0_r28.mat');
x_dotq=load('x_dot_noise0_r28.mat');
% x_dot2q=load('x_dot2.mat');
YYq=load('Y_noise0_r28.mat');
ZZq=load('Z_noise0_r28.mat');

x_ce=x_ceq.x_ce;
x_dot=x_dotq.x_dot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%数据循环
shuju_L=5000;
x_ce=x_ce(1:shuju_L,1);
x_ce_1=x_ce(1:shuju_L,1);

for i=1:750000/shuju_L
    x_ce=[x_ce;x_ce_1];
end

x_dot=x_dot(1:shuju_L,1);
x_dot_1=x_dot(1:shuju_L,1);
for i=1:750000/shuju_L
    x_dot=[x_dot;x_dot_1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_dot=x_dot2q.x_dot2;
YY=YYq.YY;
ZZ=ZZq.ZZ;
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
a=15;
b=3;
r=35;
global VAL;
y_1=1;
z_1=1;
%  [t1,H1]=mod_runge_kutta(@new_lor_lorenz,[1 1 1 0 1 0 1 0],t0,0,9120,x_ce,x_dot,x_2dot,a,b,r);
% [t1,H1]=mod_runge_kutta(@mod_lor_fun,[0 0 1 1 1 1 1 0 0 0 0.1],t0,0,120,x_ce,x_dot);
ini=[];
for i=1:10
    ini=[ini 0];
end
ini=[ini 1 1];
for i=1:5
    ini=[ini 0 1];
end

% R: 这里开始执行求解
[t1,H1]=mod_runge_kutta(@mod_lor_fun,ini,t0,0,3820,x_ce,x_dot,a,b,r,alfa,learn_rate);

figure
   plot(H1(14,:),'LineWidth',2);  % 这是预测的a
   hold on
   plot(H1(16,:),'LineWidth',2);  % 这是预测的r
   hold on
   plot(H1(18,:),'LineWidth',2);  % 这是预测的b
   xlabel('Step');
ylabel('Parameter value');
H1(14,end)
H1(16,end)
H1(18,end)

stop


% title('Lorenz D=0.1, \alpha=3, \gamma=0.01 a_{1}=15 b_{1}=3 r_{1}=35,Largest Lyapunov Exponents=1',"FontName","Times New Roman","FontSize",35,'FontWeight','bold');
title({'Reconstructed parameters using 5000 data points';'Lorenz D=0, \alpha=3, \gamma=0.01 a=10 b=8/3 r=28'},'FontName','Times New Roman','FontSize',35,'FontWeight','bold');
set(gca,'fontsize',20,'fontname','Times New Roman','FontWeight','bold');
set(gca,'FontSize',20); 
set(gca,'linewidth',3)%加粗坐标轴
xlim([0 200000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,1,1);
%  plot(H2(11,90000:93000),'LineWidth',3);
 plot(H1(11,:),'LineWidth',5);
 set(gca,'FontSize',20);
xlabel('Step','fontsize',36);
ylabel('y_{value}','fontsize',36); 
set(gca,'linewidth',5)%鍔犵矖鍧愭爣杞?
set(gca,'Xcolor',[0 0 0]);
set(gca,'Ycolor',[0 0 0]);
 hold on
%  plot(x_ce(9,90000:93000),'LineWidth',1);
 plot(YY,'linestyle','--','LineWidth',3);
legend('Reconstructed time series','True time series','FontName','Times New Roman','FontSize',30,'FontWeight','bold');
 set(gca,'XLim',[400000 403000]);
 set(gca,'YLim',[-50 50]);
%  title('Lorenz D=1',"FontName","Times New Roman","FontSize",35,'FontWeight','bold');
 subplot(2,1,2);
%   plot(H2(12,90000:93000),'LineWidth',3);
 plot(H1(12,:),'LineWidth',5);
 set(gca,'FontSize',20);
 xlabel('Step','fontsize',36);
 ylabel('z_{value}','fontsize',36); 
 set(gca,'linewidth',5)%鍔犵矖鍧愭爣杞?
 set(gca,'Xcolor',[0 0 0]);
 set(gca,'Ycolor',[0 0 0]);
 hold on
%  plot(x_ce(17,90000:93000),'LineWidth',1);
  plot(ZZ,'linestyle','--','LineWidth',3);
 set(gca,'XLim',[400000 403000]);
% title('Lorenz D=0.1',"FontName","Times New Roman","FontSize",35,'FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   figure
 plot(H1(4,:),'LineWidth',2);
 hold on
 plot(H1(3,:),'LineWidth',2);
 hold on
 plot(H1(5,:),'LineWidth',2);
  set(gca,'FontSize',20);
xlabel('Step','fontsize',36);
ylabel('value','fontsize',36); 
set(gca,'linewidth',5)%閸旂姷鐭栭崸鎰垼鏉?
set(gca,'Xcolor',[0 0 0]);
set(gca,'Ycolor',[0 0 0]);
 figure
 subplot(1,2,1);
 plot(H1(1,3000:5000),'LineWidth',3);
   set(gca,'FontSize',20);
xlabel('Step','fontsize',36);
ylabel('y_{value}','fontsize',36); 
set(gca,'linewidth',5)%閸旂姷鐭栭崸鎰垼鏉?
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
set(gca,'linewidth',5)%閸旂姷鐭栭崸鎰垼鏉?
set(gca,'Xcolor',[0 0 0]);
set(gca,'Ycolor',[0 0 0]);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
plot(g(1:450000,1),'b','LineWidth',2)
set(gca,'FontSize',20);
xlabel('Step','fontsize',36);
ylabel('absolute error','fontsize',36); 
set(gca,'linewidth',5)%閸旂姷鐭栭崸鎰垼鏉?
set(gca,'Xcolor',[0 0 0]);
set(gca,'Ycolor',[0 0 0]);
title('Error');
ylim([-300 300])

mean(H1(14,400000:end))
std(H1(14,400000:end))

mean(H1(16,400000:end))
std(H1(16,400000:end))
mean(H1(18,400000:end))
std(H1(18,400000:end))
