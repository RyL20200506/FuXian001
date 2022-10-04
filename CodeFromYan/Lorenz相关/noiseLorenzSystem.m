% Taylor cacluator nosie Lorenz System
close all;clear;clc;  % R: ��������
% for D=1:3  


for zkkk=1:2  % ?��: Ӧ����ѭ���Ĵ��� 

    % C.1 ��������
    delta=10;b=8/3;r=28;  % Lorenz��������
    N=10000;
    x0 = [-8;7;27];  % ��ʼֵ
    X_n = [x0];  % R: ����
    t0 = 0.01;  % ��λ����ʱ��
    D = 500;  % ��������
    m = 1;  % R: Ŷ���������ʱ�䴰��֮���

    % C.2. ����΢�ַ��̵Ĺ켣 ����������
    b1=sqrt(2*D*t0); % ??? ����ɶ, �����������й�
    b2=1/2*b1*t0;
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
        X_t=X_n(:,i)+[fx;fy;fz]*t0+1/2*t0*t0*df1+df2.*s2'+s1';  % R: �����һ������ȥ����״̬, ������״̬����ʱ����������.
        X_n=[X_n,X_t];  % �ϲ��µ���, ���õ�һ��3��, 10001�еľ���. �����е�����ʱ������.
    end

    % C.3 ��������õ��Ĺ켣ͼ
    A=X_n(:,1:end-1);  % ȥ�������һ��, ��Ϊ���һ��ֻ��Ϊ����ж�, ����Ϊ��֪
    subplot(1,2,1)  % 
    plot(A(1,:),A(3,:),'linewidth',1);
    xlabel('x');
    ylabel('z');
    set(gca,'FontSize',30) ;
    subplot(1,2,2)
    plot(A(1,:),'linewidth',1);
    xlabel('step');
    ylabel('x');
    set(gca,'FontSize',30) ;

    % C.4 ��ȡ�����ʹ�õĻ�����
    BASE = [];  % ? �о���仰����ɾ��
    base = baseEuqation(A');  % A������������ʱ������, A'�޸�Ϊ����������ʱ������.
    BASE = [base*t0;base*t0*t0];  % R: ����t0��Ӧԭ�ĵ�t_{unit}�ĸ���, ��ʵ����g_i��h_i, ���ǵĻ�������ͬ, ���ǵ�
    
    % C.5 (Yan)m=1,�ӵ�һ��Ԫ������,��ʼ����Ϊ1��L����2��L+1��ϵĽ��
    DBASE=[base;2*t0*base];  % R ����һ����ʽ�����˻�����, �õ���ԭ��������һ�׻�����
    BASE=BASE(:,1:N-1);  % ȥ�����һ�л����� ???����ΪʲôҪ������, ��ΪBASE�ǻ���base, base�ǻ���A, A�Ѿ���X_n��ȥ�����һ����.
    BASE_pinv=pinv(BASE);  % ׼���û�������α��

    % ??? Ϊʲô: �Ǹ���: ���=ϵ��*BASE, ����ϵ��=���/BASE��?, ��diff_x������ϵ���ĺ�����
    diff_x=(A(1,2:end)-A(1,1:end-1))*BASE_pinv;  % R ����, �ٳ��Ի�����, �����·��ܵõ�, һ������ʹ��:�������*������=�����ڵĲ��.
    diff_y=(A(2,2:end)-A(2,1:end-1))*BASE_pinv;  % 
    diff_z=(A(3,2:end)-A(3,1:end-1))*BASE_pinv;  % 
    
    % ??? DBASE��������ʲô, DX��ʾ���ǵ����ĺ�����? Ϊʲô�������������㵼��?
    DX=diff_x*DBASE;  % ��:��BASE��DBASE�����ֲ�ͬ�Ļ�����.
    DY=diff_y*DBASE;  % ���������֮��Ż��õ�, ��Ҳ��֪����û����
    DZ=diff_z*DBASE;

    % ??? ���ѵ��ǵõ�����һʱ�̵�״̬����, ����diff_xȷʵ��ʾϵ���ĺ���, ��Ϊdiff_x*BASE�͵õ��˲��.
    % ??? ����Ҳû��Ԥ��ѽ, ��ΪA(1,2:end)����������֪�Ĳ�����?
    X=A(1,1:end-1)+diff_x*BASE;  % R �õ���Ϻ�x����������
    Y=A(2,1:end-1)+diff_y*BASE;  % 
    Z=A(3,1:end-1)+diff_z*BASE;  % 

    % X=diff_x*BASE
    % Y=diff_y*BASE;
    % Z=diff_z*BASE;

    % plot(X_n(1,:),X_n(2,:))
    
    
    for i=1:N  % ??? AӦ��ֻ��3��, ������ǵ�2�м���1��, �൱��y-x, ����ÿһ��ʱ��y-x��ֵ.
        base_a(i,:)=(A(2,i)-A(1,i));  % ??? �ⲻ��y-x��, ��a��ʲô��ϵ��? base_a��ô�����?
    end
    for i=1:N
        base_b(i,:)=A(3,i);  % ??? ��ȫ�ճ�z
    end
    for i=1:N
        base_r(i,:)=A(1,i);  % ??? ��ȫ�ճ�x, ��Ϊʲô��_r, x��r��ʲô��ϵ��?
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
    
    base_ap=pinv(base_a);  % ??? ��y-x����?
    base_bp=pinv(base_b);  % 
    base_rp=pinv(base_r);
    
    for i=1:N
        cha_x(i,:)=DX(1,i);  % Why?, ΪʲôͻȻ��DXȥ��
    end
    for i=1:N
        cha_z(i,:)=A(1,i)*A(2,i)-DZ(1,i);  % ΪʲôҪ����д
    end
      for i=1:N
            cha_y(i,:)=DY(1,i)+A(2,i)+A(1,i)*A(3,i);  % ����������, �������Ŀ�������
      end
      a_value(zkkk,1)=base_ap*cha_x;  % a_value��¼������ϵ�a����ʷ.
      b_value(zkkk,1)=base_bp*cha_z;  % b_value
      r_value(zkkk,1)=base_rp*cha_y;
end

a_value(end)
b_value(end)
r_value(end)

stop
%     %������11ʽ��ȷ����
%     for i=1:size(A,2)
%         base_a(i,:)=(A(2,i)-A(1,i))*t0;  % C: base_a���Ƕ�A���������ڵ�ʱ����������;
%     end
%     base_ap=pinv(base_a); 
%     BASE_11=[base*t0; base*t0*t0];  % BASE_11��BASE����, ֻ��BASE����һ��ά�ȶ���;
% 
%     % a_value��¼Ԥ��ֵ, ��: ��Ԥ��ĵ�������ʽ��ʲô, ��ô�������е����۶�Ӧ��?
%     for i=1:size(A,2)
%        DERA(i,:)=0.5*t0*t0*(-a_value(zkkk,1)*a_value(zkkk,1)*(A(2,i)-A(1,i))+a_value(zkkk,1)*(r_value(zkkk,1)*A(1,i)-A(2,i)-A(1,i)*A(3,i)));
%     end
%     a_value(zkkk,1)=base_ap*(BASE_11'*diff_x'-DERA);  % C: �����¼�µ�aֵ
% 
%     % ����r
%     for i=1:size(A,2)
%         base_y(i,:)=A(1,i)*t0;  % Ϊɶ?
%     end
%     base_yp=pinv(base_y);
%     for i=1:size(A,2)
%        DERA_y(i,:)=0.5*t0*t0*( a_value(zkkk,1)*(A(2,i)-A(1,i))*(r_value(zkkk,1)-A(3,i))-( r_value(zkkk,1)*A(1,i)-A(2,i)-A(1,i)*A(3,i))-A(1,i)*(A(2,i)*A(1,i)-b_value(zkkk,1)*A(3,i)));
%     end
%     r_value(zkkk,1)=base_yp*(BASE_11'*diff_y'-DERA_y+t0*(A(2,:)+A(1,:).*A(3,:))');  % C: �����¼�µ�rֵ
% 
%     % ����b
%     for i=1:size(A,2)
%         base_z(i,:)=A(3,i)*t0;
%     end
%     base_zp=pinv(base_z);
%     for i=1:size(A,2)
%        DERA_z(i,:)=0.5*t0*t0*(a_value(zkkk,1)*(A(2,i)-A(1,i))*A(2,i)+(r_value(zkkk,1)*A(1,i)-A(2,i)-A(1,i)*A(3,i))*A(1,i)-b_value(zkkk,1)*(A(2,i)*A(1,i)-b_value(zkkk,1)*A(3,i)));
%     end
%     b_value(zkkk,1)=-base_zp*(BASE_11'*diff_z'-DERA_z-t0*(A(1,:).*A(2,:))');
% 
%     % �������
%     err_a(zkkk,1)=abs(a_value(zkkk,1)-delta)/delta;  % aҲ������һ�����ֽ�delta
%     err_b(zkkk,1)=abs(b_value(zkkk,1)-b)/b; 
%     err_r(zkkk,1)=abs(r_value(zkkk,1)-r)/r;  % �����µ�(λ�ڵ�zkkk�е�)r_value����ʵ��rȥ��������
%     noise_DD=0;
% 
%     for i=1:N-1
%        noise_D=(X(1,i)-A(1,i+1))^2;
%        noise_DD= noise_DD+ noise_D;
%     end
% 
%     noise_DD= noise_DD/(N-1);  % C 
%     D_D(zkkk,1)=noise_DD/(2*t0);  % C ������������
%     err_D(zkkk,1)=abs(D_D(zkkk,1)-D)/D;  % C 
    % delta=a_value;b=b_value;r=-r_value;
% end
% 
% a_err=(mean(a_value)-delta)/delta;
% b_err=(mean(b_value)-b)/b;
% r_err=(mean(r_value)-r)/r;
% D_err=(mean(D_D)-D)/D;
% 
% a_err=mean(err_a)  % C: ���������ȡ��ֵ, ��Ӧ�������Ľ����
% b_err=mean(err_b)  
% r_err=mean(err_r)
% D_err=mean(err_D)  

%  N=10000;
% % x0=[-1;3;4];
% X_n=[];
% x0=[-8;7;27];
% 
% X_n=[x0];
% t0=0.01;%��λ����ʱ��
% % D=1;%��������
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
