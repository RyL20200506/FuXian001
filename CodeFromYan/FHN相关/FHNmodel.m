clear;clc;close all;   
%% 重要参数选取
summ=zeros(4,5);
for ppr=1:20
N=2;T1=1000;t0=0.1; G=10000;   
delta=0.4;c1=0;c2=0.7;b=0.8;
initial=2*rand(1,2*N);
% initial=[-1,-1,-1,-1];
%% 产生一个N*N的矩阵，取值范围[0,1）
% W2=W;
% W1=rand(N);
% W=zeros(N);
W=rand_net(N,0.2);  % 生成随机的连接矩阵
% bbk=load('intra_A.mat');
% W=bbk.a;
for i=1:N
   for j=1:N
     if i<=j
      
         W(i,j)=0;
     end
   end
end
W=0.5*W;
% 
% W(2,1)=1;
% W(3,2)=1;
% W(4,2)=1;
% W(3,4)=0.2;
% W(5,4)=0.2;
% W(4,1)=0;
% W(4,3)=5;
% W(5,4)=0;
% W(5,2)=1;
% for i=1:N 
%     temp=0;
%     for j=1:N
%         temp=temp+W(i,j);
%     end
%     W(i,i)=-temp;
% end
% W2=roundn(W,-2);
%  save('WG10000t0.07.mat','W');  
% 绘图
%  IDS={'1','2','3','4','5','6','7','8','9','10'};
%  bg=biograph(W2,IDS);
%  set(bg.nodes,'shape','circle','color',[1,1,1],'lineColor',[0,0,0]); 
%  set(bg,'layoutType','radial'); 
%  bg.showWeights='on';
%  set(bg.nodes,'textColor',[0,0,0],'lineWidth',2,'fontsize',9);
%  set(bg,'arrowSize',2,'edgeFontSize',9);
%  get(bg.nodes,'position');
%  view(bg);
%  help biograph

% fstr=diffEquation(N);     %根据邻接矩阵生成微分方程字符串形式
% f=eval(fstr);
% [~,V1]=rk_4(f,[0,T1,t0],initial);
% [~,V2]=rk_4(f,[T1,T2,t0],V1(end,:));
% figure
% plot(V1(:,1),'b');
% xlim([0,G]);
% str=['t0=' num2str(t0) ',G=' num2str(G) 'V1.mat'];
% save(str,'V1');
%% 改进的龙格-库塔方法给系统中加噪声
fstr1=diffEquation1(N);     %根据邻接矩阵生成微分方程字符串形式
F=eval(fstr1);
Q1=0.005; Q2=0.01;
% Q1=0; Q2=0;
w1=randn(G,N);w2=randn(G,N);
X_n=[initial];
for i=1:G
     F1=F(X_n(i,:))';
     F2=F(X_n(i,:)+t0*F1(1,:) + [sqrt(2*Q1*t0)*w1(i,:),sqrt(2*Q2*t0)*w2(i,:)])';  % 后面这一项是生成白噪声？
     X1=X_n(i,:)+1/2*(F1+F2)*t0 + [sqrt(2*Q1*t0)*w1(i,:),sqrt(2*Q2*t0)*w2(i,:)];
     X_n=[X_n;X1];
end
figure(1)
plot(X_n(:,1));
set(gca,'fontsize',30,'fontname','Times New Roman');
xlim([0,G]);
% xlabel('t');
ylabel('v_1');
figure(2)
plot(X_n(:,2));
set(gca,'fontsize',30,'fontname','Times New Roman');
xlim([0,G]);
% xlabel('t');
ylabel('v_2');
% str1=['noiseFHNmodelQ1=' num2str(Q1) ',Q2=' num2str(Q2) 't0=' num2str(t0) ',G=' num2str(G) '.mat'];
X_n_s=[];
for i=1:N
    X_n_s=[ X_n_s,X_n(:,i)];
end
pmi=[];
% save(str1,'X_n_s');
% str2=['W.mat'];
% save(str2,'W');
% stop

%%%%%mutual information%%%%%
for i=1:1000
X_NN(i,:)=X_n_s(10*i,:);
end
%
for ii=1:N
    for jj=1:N
        if ii~=jj
%         for kk=1:N
        for tau=1:100
              X1=[];X2=[];X3=[];
            for i=tau:1000
                X1(i-tau+1,:)=X_NN(i,ii);
%     X3(i-tau+1,:)=X_NN(i,2);
            end
            for i=1:1000-tau+1
                X2(i,:)=X_NN(i,jj);
%     X3(i,:)=X_NN(i,1);
            end
mi(tau,:)=VectorMI(X1,X2,20);
% mi_1(tau,:)=transfer_entropy(X1,X2,X3,20);
% mi_2(tau,:)=pVectorMI(X1,X2,X3,20);
         for pp=1:N
                if pp==ii || pp==jj
                else
                    X3=[X3,X_NN(:,pp)];
                end
         end
         X3=X3(1:1000-tau+1,:);
         YYY=[X1,X2,X3];
         mi_3(tau,:)=0.2+kp_entropy(YYY,10);
% mi_4(tau,:)=kp_entropy1(X1,X2,X3,10);
        end
% for bk=2:100
%     ppmi(bk-1,1)=mi_3(1,1)-mi_3(bk,1);
% end
mean_v=mean(mi_3);
mse=0;
for i=1:size(mi_3)
    mse=mse+(mi_3(i,1)-mean_v)^2;
end
thta=(mse/(1000-1))^0.5;
MI(ii,jj)=mi_3(1,1);
MI_m(ii,jj)=max(mi);
pmi=[pmi,mi_3];
        end
    end
    
end
% end
for i=1:N
   for j=1:N
     if i<=j
      
         MI(i,j)=0;
     end
     if MI(i,j)<0.2
         MI(i,j)=0;
     end
   end
end
for tau=1:100
    X=[];Y=[];
for i=tau:1000
   X2(i-tau+1,:)=X_NN(i,4);
   X3(i-tau+1,:)=X_NN(i,5);
end
for i=1:1000-tau+1
     X1(i,:)=X_NN(i,2);
%      X3(i,:)=X_NN(i,2);
end
mi_3(tau,:)=pVectorMI(X2,X1,X3,20);
end
T_AU=[-99:100];
mi_4=[];
for i=100:-1:1
   mi_4=[mi_4;mi_3(i,1)];
end
MI=[mi_4;mi_2];
T_TAU=[-99:100]';
plot(T_TAU,MI);
set(gca,'FontSize',30) ;
xlabel('tau');
ylabel('PMI');
axis([-inf inf,0,0.2]);
%%%%%%
num=size(X_n_s);
NUM=num(1,1);
node_num=N;
b_value=zeros(NUM,node_num);
for j=1:node_num
VV_value=double(X_n_s(:,j));
mean_v=mean(VV_value);
mse=0;

for i=1:NUM
    mse=mse+(VV_value(i,1)-mean_v)^2;
end

thta=(mse/(NUM-1))^0.5;
thta=2*thta;

for i=1:NUM
    if VV_value(i,1)>thta
        b_value(i,j)=1;
   
    end
end
end
% for j=1:node_num
% t_value=[];
% for i=1:NUM
%     if b_value(i,j)==1
%        t_value=[t_value;i];
%     end
% end
% 
% fid = fopen(['C:\Users\yan\Desktop\非线性拟合——王娟娟\fhn_',num2str(j),'.txt'],'wt'); 
% fprintf(fid,'%10f\n',t_value); 
% fclose(fid);
% end
CXY=[];
W_TAU=100;
for s=1:node_num
    for u=s:node_num
        if s~=u
        tC_xy=[];
        C_xy=[];
        s_sum=sum(b_value(:,s));
        u_sum=sum(b_value(:,u));
for tau=-W_TAU:0
    k_mean=0;
for i=-tau+1:NUM
    k_mean=k_mean+b_value(i,s)*b_value(i+tau,u);
end
tC_xy=[tC_xy;tau];
C_xy=[C_xy;(1/(s_sum*u_sum)^0.5)*k_mean];
% C_xy(tau,2)=(1/(sum(b_value(:,1))*sum(b_value(:,2)))^0.5)*k_mean;
end

for tau=1:W_TAU
    k_mean=0;
for i=1:NUM-tau
    k_mean=k_mean+b_value(i,s)*b_value(i+tau,u);
end
tC_xy=[tC_xy;tau];
C_xy=[C_xy;(1/(s_sum*u_sum)^0.5)*k_mean];

end

C_xy=C_xy-mean(C_xy);
CXY=[CXY,C_xy];
if max(C_xy)>abs(min(C_xy))
    cc_out(s,u)=max(C_xy);
    [m,p]=max(C_xy);
%     if p<70
%         cc_out(s,u)=0;
%     end
    if p>=1000
        cc_out_fangxaing(s,u)=1;
    else
        cc_out_fangxaing(s,u)=-1;
    end
else 
   cc_out(s,u)=min(C_xy);
       [m,p]=min(C_xy);
    if p>=1000
        cc_out_fangxaing(s,u)=1;
    else
        cc_out_fangxaing(s,u)=-1;
    end
   
end
    
        end
    end
end
summ=summ+cc_out;
% plot(CXY(:,1))
% hold on
% plot(CXY(:,2))
% hold on
% plot(CXY(:,3))
% hold on
% plot(CXY(:,4))
% hold on
% plot(CXY(:,5))
% hold on
% plot(CXY(:,6))
end
summ=summ/20;
for i=1:4
    for j=1:5
       if summ(i,j)<0.2
           summ(i,j)=0;
       end
    end
end