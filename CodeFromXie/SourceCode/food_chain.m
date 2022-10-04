function [dR,dC,dP]=food_chain(R,C,P,x_c,y_c,z_c)
dR=x_c*(C-R);
dC=z_c*R-C-R*P;
dP=R*C-y_c*P;
end
% function [dR,dC]=food_chain(t,R,C,b1,k,M1,G,f) % define the function
% % global t
% dR=C;
% dC=-b1*C-k*R+M1*sin(R)+G*cos(2*pi*f*t);
% end
% function [dR,dC]=food_chain(t,R,C,b1,b2,k,M1,M2,G,f) % define the function
% % global t
% dR=C;
% dC=-b1*C-b2*C.^2-k*R+M1*sin(R)+M2*sin(2*R)+G*cos(2*pi*f*t);
% end