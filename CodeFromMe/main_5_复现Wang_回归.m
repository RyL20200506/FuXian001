% ����: Wang - �ع�����ѧ: �����Իع������

close all;clear;clc;

% ������֪����
y0 = [-8,7,27];  % ��ʼֵ
[t, real_trajectory] = ode45('fun_0_lorenz', [0,100], y0);  % ����֮ǰ�ĺ��������ʵ�켣

% һЩ����
t_unit = 0.01; 

% ׼�����ݼ�
% (1)ȡ������
X = real_trajectory(:,1);  % ȡ��������, ʱ������������
Y = real_trajectory(:,2);
Z = real_trajectory(:,3);

% (2)׼��: ������Ҫ�Ļ�����
XY = X.*Y; 
XZ = X.*Z; 
YZ = Y.*Z; 
XYZ = XY.*Z;

% (3)�ϲ�
about_all_x = [ones(size(X)), X, X*t_unit, Y*t_unit, Z*t_unit, XY*t_unit, XZ*t_unit, YZ*t_unit, XYZ*t_unit, X*(t_unit^2), Y*(t_unit^2), Z*(t_unit^2), XY*(t_unit^2), XZ*(t_unit^2), YZ*(t_unit^2), XYZ*(t_unit^2)];
about_all_y = [ones(size(Y)), Y, X*t_unit, Y*t_unit, Z*t_unit, XY*t_unit, XZ*t_unit, YZ*t_unit, XYZ*t_unit, X*(t_unit^2), Y*(t_unit^2), Z*(t_unit^2), XY*(t_unit^2), XZ*(t_unit^2), YZ*(t_unit^2), XYZ*(t_unit^2)];
about_all_z = [ones(size(Z)), Z, X*t_unit, Y*t_unit, Z*t_unit, XY*t_unit, XZ*t_unit, YZ*t_unit, XYZ*t_unit, X*(t_unit^2), Y*(t_unit^2), Z*(t_unit^2), XY*(t_unit^2), XZ*(t_unit^2), YZ*(t_unit^2), XYZ*(t_unit^2)];

% ѵ��
known_all_x = about_all_x(1:end-1,:);  % ��Ϊ���һ��ֻ�䵱label
need_all_x = X(2:end);  % ��Ϊ��һ��ֻ�䵱feature
para_x = regress(need_all_x, known_all_x)  % ���Իع�

known_all_y = about_all_y(1:end-1,:);  % ��Ϊ���һ��ֻ�䵱label
need_all_y = Y(2:end);  % ��Ϊ��һ��ֻ�䵱feature
para_y = regress(need_all_y, known_all_y)  % ���Իع�

known_all_z = about_all_z(1:end-1,:);  % ��Ϊ���һ��ֻ�䵱label
need_all_z = Z(2:end);  % ��Ϊ��һ��ֻ�䵱feature
para_z = regress(need_all_z, known_all_z)  % ���Իع�

% -----------------------------------------

% �鿴ѵ��Ч��
% ��������X�Ľ��
the_known_x = known_all_x(1,:);
the_knwon_y = known_all_y(1,:);
the_known_z = known_all_z(1,:);

% 
X_predict_history = [];  % ���� ��¼���
Y_predict_history = []; 
Z_predict_history = []; 

for i=1:size(known_all_x)
    x = the_known_x*para_x;  % 
    y = the_knwon_y*para_y;
    z = the_known_z*para_z;
    
    X_predict_history = [X_predict_history; x];  % �кϲ� 
    Y_predict_history = [Y_predict_history; y]; 
    Z_predict_history = [Z_predict_history; z]; 
    
    % ���»�����
    the_known_x = [1, x, x*t_unit, y*t_unit, z*t_unit, (x*y)*t_unit, (x*z)*t_unit, (y*z)*t_unit, (x*y*z)*t_unit, x*(t_unit^2), y*(t_unit^2), z*(t_unit^2), x*y*(t_unit^2), x*z*(t_unit^2), y*z*(t_unit^2), x*y*z*(t_unit^2)];
    the_known_y = [1, y, x*t_unit, y*t_unit, z*t_unit, (x*y)*t_unit, (x*z)*t_unit, (y*z)*t_unit, (x*y*z)*t_unit, x*(t_unit^2), y*(t_unit^2), z*(t_unit^2), x*y*(t_unit^2), x*z*(t_unit^2), y*z*(t_unit^2), x*y*z*(t_unit^2)];
    the_known_z = [1, z, x*t_unit, y*t_unit, z*t_unit, (x*y)*t_unit, (x*z)*t_unit, (y*z)*t_unit, (x*y*z)*t_unit, x*(t_unit^2), y*(t_unit^2), z*(t_unit^2), x*y*(t_unit^2), x*z*(t_unit^2), y*z*(t_unit^2), x*y*z*(t_unit^2)];
end

% ��ͼ
plot(X_predict_history, 'r');
hold on

% y0 = [1,1,1];  % ��ʼֵ
% [t, real_trajectory] = ode45('fun_0_lorenz', [0,200], y0);  % ����֮ǰ�ĺ��������ʵ�켣
% plot(real_trajectory(:,1), '*');
% hold on
% 
% % -----------------------------------------------
% % ò�ƺܺ�, ��ʵ�Ǵ���Ľ��
% the_known_x = known_all_x(1,:);
% the_knwon_y = known_all_y(1,:);
% the_known_z = known_all_z(1,:);
% 
% % 
% X_predict_history = [];  % ���� ��¼���
% Y_predict_history = []; 
% Z_predict_history = []; 
% 
% for i=1:size(known_all_x)
%     x = the_known_x*para_x;  % 
%     y = the_knwon_y*para_y;  % 
%     z = the_known_z*para_z;  % 
%     
%     X_predict_history = [X_predict_history; x];  % �кϲ� 
%     Y_predict_history = [Y_predict_history; y];  % 
%     Z_predict_history = [Z_predict_history; z];  
%     
%     % ���»�����
%     the_known_x = [1, x, x*t_unit, y*t_unit, z*t_unit, (x*y)*t_unit, (x*z)*t_unit, (y*z)*t_unit, (x*y*z)*t_unit, x*(t_unit^2), y*(t_unit^2), z*(t_unit^2), x*y*(t_unit^2), x*z*(t_unit^2), y*z*(t_unit^2), x*y*z*(t_unit^2)];
%     the_known_y = [1, y, x*t_unit, y*t_unit, z*t_unit, (x*y)*t_unit, (x*z)*t_unit, (y*z)*t_unit, (x*y*z)*t_unit, x*(t_unit^2), y*(t_unit^2), z*(t_unit^2), x*y*(t_unit^2), x*z*(t_unit^2), y*z*(t_unit^2), x*y*z*(t_unit^2)];
%     the_known_z = [1, z, x*t_unit, y*t_unit, z*t_unit, (x*y)*t_unit, (x*z)*t_unit, (y*z)*t_unit, (x*y*z)*t_unit, x*(t_unit^2), y*(t_unit^2), z*(t_unit^2), x*y*(t_unit^2), x*z*(t_unit^2), y*z*(t_unit^2), x*y*z*(t_unit^2)];
% end
% 
% % ��ͼ
% plot(X_predict_history, 'r');
% hold on
% 
% 
