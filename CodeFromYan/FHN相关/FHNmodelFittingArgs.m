clear;clc;close all; 
%noiseFHNmodelQ1=0.005,Q2=0.01t0=0.01,G=10000 误差很大
%noiseFHNmodelQ1=0.005,Q2=0.01t0=0.02,G=10000 误差很大
%noiseFHNmodelQ1=0.005,Q2=0.01t0=0.05,G=10000 误差很大
%t=0.1时的误差最小，但随着噪声的增大误差也逐渐增大
%% 重要参数选取
N=10;
T1=800;t0=0.08;G=T1/t0;
A1=load('noiseFHNmodelQ1=0.005,Q2=0.05t0=0.08,G=10000.mat');
B=A1.X_n;
V=[];W=[];RES_W=[];RES_V=[];DV=[];DW=[];temp_c=[];WW=[];C2=[];res_c1=[];res_c2=[];
res_b=[];res_delta=[];
for i=1:N
    base=baseEuqationV(B,N,i);  
  
    Base_v=[[base;B(:,1:10)']*t0;base*t0*t0];   
    base_w=baseEuqationV(B,N,i);
    Base_w=[base_w*t0;base_w*t0*t0];

%     Base_v_T=Base_v*Base_v';

%     Base_v_pinv=pinv(Base_v_T);
    
%     ceshi =Base_v_pinv*Base_v_T;
    Base_v_pinv=pinv(Base_v(:,1:end-1));
    Base_w_pinv=pinv(Base_w(:,1:end-1));

%     res_v=Base_v_pinv*Base_v(:,1:end-1)*(B(2:end,i)-B(1:end-1,i));
    res_v=(B(2:end,i)-B(1:end-1,i))'*Base_v_pinv;%计算出模型分量V的假设系数
    res_w=(B(2:end,10+i)-B(1:end-1,10+i))'*Base_w_pinv;%计算出模型分量W的假设系数
    RES_V=[RES_V;res_v]; RES_W=[RES_W;res_w];
    v1=B(1:end-1,i)'+res_v*Base_v(:,1:end-1);
    w1=B(1:end-1,10+i)'+res_w*Base_w(:,1:end-1);
    V=[V;v1];W=[W;w1];
    Dbase_v=[base;B(:,1:10)';2*base*t0];
    Dbase_w=[base_w;2*base_w*t0];
    Dv=res_v*Dbase_v(:,1:end-1);%根据得到的系数得到dot v
    Dw=res_w*Dbase_w(:,1:end-1);
    DV=[DV;Dv]; DW=[DW;Dw]; 
    NK=size(base);
   
%     vv=res_v(NK(1,1)+1:NK(1,1)+N,1); % W_ij拉普拉斯矩阵第i行的数据，表示其他节点对第i个节点的影响系数
     vv=res_v(1,NK(1,1)+1:NK(1,1)+N);
%       vv=res_v(1,10:19);
%     vv(i,1)=vv(i,1)-1;  %要提出每一个vi放在前面
      vv(1,i)=vv(1,i)-1; 
      WW=[WW;vv];%要求解的拉普拉斯矩阵W_ij
     %求解fhn模型中的系数
      temp_c1=res_v(1,1);res_c1=[res_c1,temp_c1];
      temp_delta=res_w(1,2);res_delta=[res_delta,temp_delta];
      temp_b=-(res_w(1,3))/temp_delta;res_b=[res_b,temp_b];
      temp_c2=(res_w(1,1))/temp_delta;res_c2=[res_c2,temp_c2];
end
% plot(V(1,:))
% hold on
% plot(B(:,1))
c1=mean(res_c1); 
delta=mean(res_delta);
c2=mean(res_c2);
b=mean(res_b);
DW_M=mean(DW);V_M=mean(V);W_M=mean(W);
temp_c2=[c2];temp_b=[b];temp_delta=[];fin_res=[];
for i=1:50
    t1_delta=DW_M/(V_M-temp_b(1,i)*W_M+temp_c2(1,i));
    temp_delta=[temp_delta,t1_delta];
    t1_c2=mean(DW_M/t1_delta+temp_b(1,i)*W_M-V_M);
    t1_b=(V_M-DW_M/t1_delta+t1_c2)/W_M;
    temp_b=[temp_b,t1_b];temp_c2=[temp_c2,t1_c2];
end
t1_delta=DW_M/(V_M-temp_b(1,end)*W_M+temp_c2(1,end));
temp_delta=[temp_delta,t1_delta];
fin_res=[temp_delta;temp_c2;temp_b];
COFF=[fin_res(:,end);c1];

err_delta=abs(COFF(1,1)-0.08)/COFF(1,1)
err_c2=abs(COFF(2,1)-0.7)/COFF(2,1)
C1_yu=c1;
for repeat=1:100
    
for j=1:N
%     for i=1:10000
%         mean_yuv(j,i)=B(i,j)+RES_V(j,1:20)*Base_v(1:20,i)*t0+RES_V(j,21:30)*Base_v(21:30,i)*t0*t0;
%     end
    for i=1:10000
        
        mean_zhenv(j,i)=B(i,j)+(B(i,j)-1/3*power(B(i,j),3)-B(i,j+N))*t0;
        w_s=0;
        if repeat==1
            for zk=1:N
                w_s=w_s+WW(j,zk)*B(i,zk);
            end
        else
            for zk=1:N
                w_s=w_s+w_yu(j,zk)*B(i,zk);
            end
        end
        v_i(j,i)=B(i,j)-1/3*power(B(i,j),3)-B(i,j+N)+C1_yu+w_s;
        w_i(j,i)=COFF(1,1)*(B(i,j)-COFF(3,1)*B(i,j+N)+COFF(2,1));
        mean_second(j,i)=t0^2/2*(v_i(j,i)*(1-B(i,j)^2)+w_i(j,i)*-1);
        mean_vi(j,i)=mean_zhenv(j,i)+ mean_second(j,i);
    end
end
% mean_yuv=mean_yuv';
% mean_zhenv=mean_zhenv';

for i=1:10000
% for j=1:N

Base_vq(i,:)=[1,B(i,1:N)]*(-t0);
%        Base_vq(i,:)=B(i,1:N)*(-t0);
% Base_vq=[Base_vq(i,:),[1,B(i,1:N)]*(-t0)];

% end
end
% Base_vq_T=Base_vq'*Base_vq;
% det(Base_vq_T)
% Base_vq_p=pinv(Base_vq_T);
% ceshi =Base_vq_p*Base_vq_T;
Base_vq_p=pinv(Base_vq);
chazhi=mean_vi-V;
for i=1:N
    rev_q(:,i)=Base_vq_p*chazhi(i,:)';
end
w_yu=rev_q(2:11,:);
C1_yu=mean(rev_q,2);
C1_yu=C1_yu(1,1);
end

% % +2*t0^2*(1-B(i+1,j)^2)
% mean_d=mean(D_NOISE,2);
% mean(mean_d);
for j=1:N
% j=1;
% i=105;
    for i=1:10000
        D(j,i)= (B(i+1,j)-V(j,i))^2;
        D_NOISE(j,i)= D(j,i)/(2*t0+2*t0^2*(1+w_yu(j,j)-B(i+1,j)^2));
    end
end
Q_NOISE=mean(D_NOISE);
Q_NOISE=mean(Q_NOISE);
err_D=abs(Q_NOISE-0.0025)/Q_NOISE

% for i=1:N
%     for j=1:N
%         s=s+V(j,i)*WW(j,i);
%     end
% end


% save('FHN_COFF_Q1=0.0025,Q2=0.05t0=0.08,G=10000.mat','COFF');
% save('FHN_WW_Q1=0.0025,Q2=0.05t0=0.08,G=10000.mat','WW');
%拟合后的W_ij和原W_ij的比较
W=load('W.mat');
W=W.W;

te1=reshape(W,1,100);
te2=reshape(WW,1,100);
te3=reshape(rev_q(2:11,:)',1,100);
% plot(te2(1,:),'b');
% hold on
% plot(te1(1,:),'r');
% subplot(1,2,1);
% x=-7:1;
% % figure (2)
% plot(x,x,'b');
% hold on
% plot(te1(1,:),te2(1,:),'square');
% set(gca,'fontsize',30,'fontname','Times New Roman');
% xlim([-7,1]);ylim([-7,1]);
% xlabel('W_{ij}');
% ylabel('W^’_{ij}');
% % title('W_{ij} t0=0.1');
% subplot(1,2,2);
x=-7:1;
% figure (2)
plot(x,x,'b','LineWidth',3);
hold on
plot(te1(1,:),te3(1,:),'.','markersize',30);
% title(['推断网络连接强度W^’_{ij}和真实网络连接强度W_{ij}的对比'],'fontname','微软雅黑','FontWeight','bold');
set(gca,'fontsize',20,'fontname','Times New Roman','FontWeight','bold');
set(gca,'FontSize',20); 
set(gca,'linewidth',3)%加粗坐标轴
set(gca, 'XTick', [-7 -6 -5 -4 -3 -2 -1 0 1])
axis square
xlim([-7,1]);ylim([-7,1]);
xlabel('Real network connection strength W_{ij}','fontname','微软雅黑','FontWeight','bold');
ylabel('Infer network connection strength W^’_{ij}','fontname','微软雅黑','FontWeight','bold');
set(get(gca,'title'),'fontname','微软雅黑')
set(0,'defaultAxesFontName','<微软雅黑>');
te11=te1';
te22=te2';
save('te12','te2')
save('te11','te1')


