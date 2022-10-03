% ����Wang ������Yanʦ�ָ��Ĵ���������
clear
clc

rng(1)
% ����׼�� - ��1: ����������, ֱ���������������
% step = 0.01;  % ��ϲ���Ϊstep, ��ʵ�켣�Ĳ���Ϊstep/2  (��Ϊ:RungeKutta���õ�1/2�Ĳ���ȥ����)
% t_u = step;
% b_fit = 100;
% b_test = 20;
% real_b = b_fit+b_test;  % ��ʵ�켣�Ľ���ʱ�� 
% real_a = 0;  % ��ʵ�켣����ʼʱ�� 
% y0 = [-8,7,27];  % ��ʵ�켣�ĳ�ʼ�� 
% interp_step = step  % step �� step/2: #ע: ��������step/2, ���к����ص�����
% interp_a = real_a + interp_step;  % ��ֵF�Ŀ�ʼʱ��  ȡstep/2�����Ĳ�ֺ�, ǰ�����1��ʱ�̵�����;
% interp_b = real_b - interp_step;  % ��ֵF�Ľ���ʱ��
% [t,X_n] = fun_8_RungeKutta(@fun_0_lorenz,real_a,interp_step,real_b,y0);  % ���켣: �̶�����: ������������ò�����һ��
% X = X_n(1:3, 1:10000)';

% ����׼�� - ��2: ����������
delta=10;b=8/3;r=28;  % Lorenz��������
N_fit = 10000
N_test = 2000;
N=N_fit+N_test;
x0 = [-8;7;27];  % ��ʼֵ
X_n = [x0];  % R: ����
t_u = 0.01;  % ��λ����ʱ��
D = 0.01;  % ��������
m = 1;  % R: Ŷ���������ʱ�䴰��֮���
b1=sqrt(2*D*t_u); % ׼��ʹ�ù�ʽA3
b2=1/2*b1*t_u;
b3=1/sqrt(3)*b2;   
w1=randn(N,3);w2=randn(N,3);  % R: ���ɱ�׼������
for i=1:N 
    s1=b1*w1(i,:);  % 
    s2=b2*w1(i,:)+b3*w2(i,:);  % 
    df2=[0;(r-X_n(3,i)-1-X_n(1,i));(X_n(2,i)+X_n(1,i)-b)];  %

    % R: ������Lorenz����, Ӧ����ԭ����, ��������Ϸ���.
    fx=delta*(X_n(2,i)-X_n(1,i));  %
    fy=r*X_n(1,i)-X_n(2,i)-X_n(1,i)*X_n(3,i);  %
    fz=X_n(1,i)*X_n(2,i)-b*X_n(3,i);

    % R: ��֮Ӧ����������һ��X_n
    df1=[-delta*fx+delta*fy; (r-X_n(3,i))*fx-fy-X_n(1,i)*fz; X_n(2,i)*fx+X_n(1,i)*fy-b*fz];  % 
    X_t=X_n(:,i)+[fx;fy;fz]*t_u+1/2*t_u*t_u*df1+df2.*s2'+s1';  % R: �����һ������ȥ����״̬, ������״̬����ʱ����������.
    X_n=[X_n,X_t];  % �ϲ��µ���, ���õ�һ��3��, 10001�еľ���. �����е�����ʱ������.
end
X = X_n(:,1:N_fit)';

% ����Ԥ�����
base = fun_13_baseEuqation_my(X(1:end-1, :));  % (1)����Сbase����, �����������һ�е�����!
BASE = [base*t_u, base*t_u*t_u];  % ���ֹ��췽ʽ��û����ȫ�ø���˼·, ���������˺���Ԥ��.
BASE_pinv = pinv(BASE);  %
diff_x = BASE_pinv*(X(2:end,1)-X(1:end-1,1));  % R ����, �ٳ��Ի�����, �����·��ܵõ�, һ������ʹ��:�������*������=�����ڵĲ��.
diff_y = BASE_pinv*(X(2:end,2)-X(1:end-1,2));
diff_z = BASE_pinv*(X(2:end,3)-X(1:end-1,3));

% (ʧ��!)����Ԥ��: (1)��һ����ʼ�� (2)����Ԥ����һ����ʼ��
a_X = X(1,:)
predict_history = [a_X]
for i=1:3000
    a_base = fun_13_baseEuqation_my(a_X);
    a_BASE = [a_base*t_u, a_base*t_u*t_u];
    next_x = a_BASE * diff_x;
    next_y = a_BASE * diff_y;
    next_z = a_BASE * diff_z;
    % ��������
    a_X = [next_x, next_y, next_z]
    predict_history = [predict_history; a_X];  % ������������Կ���ʧ�ܵĽ��.
end

% �ȹ���a,b,r
DBASE = [base, 2*t_u*base];  % ԭ��������ʱ����
DX = DBASE*diff_x;  % ��һʱ�̵�dot{x}
DY = DBASE*diff_y;  
DZ = DBASE*diff_z;

% ��: ����lorenz��ʽ, ׼������֮�����
for i=1:size(DX,1)
    cha_x(i,:)=DX(i,1);  % lorenz ��a�Ĳ�����֮�����
end
for i=1:size(DX,1)
    cha_y(i,:)=DY(i,1)+X(i,2)+X(i,1)*X(i,3);  % lorenz ��r�Ĳ�����֮�����
end
for i=1:size(DX,1)
    cha_z(i,:)=X(i,1)*X(i,2)-DZ(i,1);  % lorenz ��b�Ĳ�����֮�����
end
% ��: ��ȡ����ʱ�̲���ǰ��ϵ��
for i=1:size(DX,1)  % 
    base_a(i,:)=(X(i,2)-X(i,1));  % lorenz ����a��ϵ��y-x
end
for i=1:size(DX,1)
    base_b(i,:)=X(i,3);  % lorenz ����b��ϵ��z
end
for i=1:size(DX,1)
    base_r(i,:)=X(i,1);  % lorenz ����x��ϵ��r
end
% ��: ׼��������
base_ap=pinv(base_a);
base_bp=pinv(base_b);
base_rp=pinv(base_r);
% ��: ���������� �õ������ع�ֵ
a_value(1)=base_ap*cha_x  % ��ϵ���������, �ɵ�a��ֵ
b_value(1)=base_bp*cha_z  % ͬ��
r_value(1)=base_rp*cha_y  % ͬ��


% ����Ԥ��: ��(1) �����������
% ���: ����ʵ����Ҳ������������ɵ������, �������ǰԤ��500������
% para = [a_value(1), b_value(1), r_value(1)]
% step = 0.01;  % ��ϲ���Ϊstep, ��ʵ�켣�Ĳ���Ϊstep/2  (��Ϊ:RungeKutta���õ�1/2�Ĳ���ȥ����)
% h = step;
% a = 0;
% b_fit = 0;
% b_test = 20;
% real_b = b_fit + b_test;  % ��ʵ�켣�Ľ���ʱ��
% real_a = 0;  % ��ʵ�켣����ʼʱ��
% y0 = X(end,:);  % ��ʵ�켣�ĳ�ʼ��
% n=floor((real_b-real_a)/h);       %����
% x(1)=a;                 %ʱ�����
% y=[];
% y(:,1)=y0;              %����ֵ������������������Ҫע��ά��
% for i=1:n               %�����������������ֵ���
%     x(i+1)=x(i)+h;
%     k1=fun_0_lorenz_para(x(i),y(:,i), para);
%     k2=fun_0_lorenz_para(x(i)+h/2,y(:,i)+h*k1/2, para);
%     k3=fun_0_lorenz_para(x(i)+h/2,y(:,i)+h*k2/2, para);
%     k4=fun_0_lorenz_para(x(i)+h,y(:,i)+h*k3, para);
%     y(:,i+1)=y(:,i)+h*(k1+2*k2+2*k3+k4)/6;
% end
% X_pred = y';


% ����Ԥ��: ��(2) ����a, b, r��value, ���������һ��ֵ����Ԥ��2000��
% ���: ����ܹ���ǰԤ��300������
delta=a_value(1); b=b_value(1); r=r_value(1);  % Lorenz��������
N_test = 2000;
N=N_test;
x0 = X(end,:);  % ��ʼֵ
X_pred = [x0'];  % R: ����
t_u = 0.01;  % ��λ����ʱ��
D = 0;  % ��������
m = 1;  % R: Ŷ���������ʱ�䴰��֮���
b1=sqrt(2*D*t_u); % ׼��ʹ�ù�ʽA3
b2=1/2*b1*t_u;
b3=1/sqrt(3)*b2;   
w1=randn(N,3);w2=randn(N,3);  % R: ���ɱ�׼������
for i=1:N 
    s1=b1*w1(i,:);  % 
    s2=b2*w1(i,:)+b3*w2(i,:);  % 
    df2=[0;(r-X_pred(3,i)-1-X_pred(1,i));(X_pred(2,i)+X_pred(1,i)-b)];  %

    % R: ������Lorenz����, Ӧ����ԭ����, ��������Ϸ���.
    fx=delta*(X_pred(2,i)-X_pred(1,i));  %
    fy=r*X_pred(1,i)-X_pred(2,i)-X_pred(1,i)*X_pred(3,i);  %
    fz=X_pred(1,i)*X_pred(2,i)-b*X_pred(3,i);

    % R: ��֮Ӧ����������һ��X_n
    df1=[-delta*fx+delta*fy; (r-X_pred(3,i))*fx-fy-X_pred(1,i)*fz; X_pred(2,i)*fx+X_pred(1,i)*fy-b*fz];  % 
    X_t=X_pred(:,i)+[fx;fy;fz]*t_u+1/2*t_u*t_u*df1+df2.*s2'+s1';  % R: �����һ������ȥ����״̬, ������״̬����ʱ����������.
    X_pred=[X_pred,X_t];  % �ϲ��µ���, ���õ�һ��3��, 10001�еľ���. �����е�����ʱ������.
end
X_pred = X_pred';

% ��ͼ
figure
plot(X_n(1,:)','linewidth',1);
hold on
t_range = 10000:1:12000;
plot(t_range, X_pred(:,1),'--', 'linewidth',1);
xlim([9500 11000]);


