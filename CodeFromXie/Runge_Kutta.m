function [R,C,P]=Runge_Kutta(dt,n,R,C,P,x_c,y_c,z_c) % R,C,P输入为初值
for k =1:n-1
    [k_1,k_2,k_3]=food_chain(R(k),C(k),P(k),x_c,y_c,z_c);

    [k_4,k_5,k_6]=food_chain(R(k)+dt/2*k_1,C(k)+dt/2*k_2,P(k)+dt/2*k_3,x_c,y_c,z_c);
    [k_7,k_8,k_9]=food_chain(R(k)+dt/2*k_4,C(k)+dt/2*k_5,P(k)+dt/2*k_6,x_c,y_c,z_c);
    [k_10,k_11,k_12]=food_chain(R(k)+dt*k_7,C(k)+dt*k_8,P(k)+dt*k_9,x_c,y_c,z_c);
    R(k+1)=R(k)+dt/6*(k_1+2*k_4+2*k_7+k_10);
    C(k+1)=C(k)+dt/6*(k_2+2*k_5+2*k_8+k_11);
    P(k+1)=P(k)+dt/6*(k_3+2*k_6+2*k_9+k_12);
end
end
% function [R,C,t]=Runge_Kutta(dt,n,t,R,C,b1,k,M1,G,f) % R,C,P输入为初值
% for i =1:n-1
%     [k_1,k_2]=food_chain(t(i),R(i),C(i),b1,k,M1,G,f); % f(x) and g(x)
%     [k_4,k_5]=food_chain(t(i)+dt/2,R(i)+dt/2*k_1,C(i)+dt/2*k_2,b1,k,M1,G,f);
%     [k_7,k_8]=food_chain(t(i)+dt/2,R(i)+dt/2*k_4,C(i)+dt/2*k_5,b1,k,M1,G,f);
%     [k_10,k_11]=food_chain(t(i)+dt,R(i)+dt*k_7,C(i)+dt*k_8,b1,k,M1,G,f);
%     R(i+1)=R(i)+dt/6*(k_1+2*k_4+2*k_7+k_10);
%     C(i+1)=C(i)+dt/6*(k_2+2*k_5+2*k_8+k_11);
%     t(i+1)=t(i)+dt;
% end   
% % R=R_0;
% % C=C_0;
% end  
% function [R,C,t]=Runge_Kutta(dt,n,t,R,C,b1,b2,k,M1,M2,G,f) % R,C,P输入为初值
% for i =1:n-1
%     [k_1,k_2]=food_chain(t(i),R(i),C(i),b1,b2,k,M1,M2,G,f); % f(x) and g(x)
%     [k_4,k_5]=food_chain(t(i)+dt/2,R(i)+dt/2*k_1,C(i)+dt/2*k_2,b1,b2,k,M1,M2,G,f);
%     [k_7,k_8]=food_chain(t(i)+dt/2,R(i)+dt/2*k_4,C(i)+dt/2*k_5,b1,b2,k,M1,M2,G,f);
%     [k_10,k_11]=food_chain(t(i)+dt,R(i)+dt*k_7,C(i)+dt*k_8,b1,b2,k,M1,M2,G,f);
%     R(i+1)=R(i)+dt/6*(k_1+2*k_4+2*k_7+k_10);
%     C(i+1)=C(i)+dt/6*(k_2+2*k_5+2*k_8+k_11);
%     t(i+1)=t(i)+dt;
% end   
% % R=R_0;
% % C=C_0;
% end  
% function [theta1,theta2,theta3]=Runge_Kutta(dt,n,W,N,KK,theta1,theta2,theta3) % R,C,P输入为初值
% for k =1:n-1
%     [k_1,k_2,k_3]=K(W,N,KK,theta1(k),theta2(k),theta3(k));
%     [k_4,k_5,k_6]=K(W,N,KK,theta1(k)+dt/2*k_1,theta2(k)+dt/2*k_2,theta3(k)+dt/2*k_3);
%     [k_7,k_8,k_9]=K(W,N,KK,theta1(k)+dt/2*k_4,theta2(k)+dt/2*k_5,theta3(k)+dt/2*k_6);
%     [k_10,k_11,k_12]=K(W,N,KK,theta1(k)+dt*k_7,theta2(k)+dt*k_8,theta3(k)+dt*k_9);
%     theta1(k+1)=theta1(k)+dt/6*(k_1+2*k_4+2*k_7+k_10);
%     theta2(k+1)=theta2(k)+dt/6*(k_2+2*k_5+2*k_8+k_11);
%     theta3(k+1)=theta3(k)+dt/6*(k_3+2*k_6+2*k_9+k_12);
% end
% end