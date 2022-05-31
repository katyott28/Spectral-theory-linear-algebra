function [v0,v1,v2]=Spec5(v)

% Spec5(v)
%   
%   Input a column vector from R^10 v=[a1;a2;a3;a4;a5;a6;a7;a8;a9;a10]
%   Outputs three vectors v0, v1, v2 so that 
%   v=v0+v1+v2
%   with each vi an effects vector
%
%   Example: Find 0, 1, and 2 effects vector for [2;12;11;6;17;8;4;24;20;6]
%
%   v=[2;12;11;6;17;8;4;24;20;6]
%   Spec5(v) gives v0=[11;11;11;11;11;11;11;11;11;11]; 
%   v1=[ ]; 
%   v2=[ ]
%   
%   Also produces "Pure 1 alignment" calculations.  Normalized projection
%   of v1 into M_1 space; Normalized projection of a Rows 1-4 of "counting"
%   matrix T1; dot product first step with second step results together to
%   create "alignment" values
%
%   Finally, produces "Pure 2 alignment" calculations.  Normalized
%   projection of v2 into M_2 space; Normalized projection of a Rows 1-6 of identity
%   matrix T2; dot product first step with second step results together to
%   create "alignment" values

%The Matrices for this context, with the transpose of T0 at the end
T1=[1 1 1 1 0 0 0 0 0 0; 1 0 0 0 1 1 1 0 0 0; 0 1 0 0 1 0 0 1 1 0; 0 0 1 0 0 1 0 1 0 1; 0 0 0 1 0 0 1 0 1 1];
T2=eye(10);
z0=(1/10)*ones(10,1);
%Calculation of Effects vectors using projections

u1=dot(v,z0)/dot(z0,z0)*z0;
v0=u1;

u2=T1(2,:)'-(dot(u1, T1(2,:)')/dot(u1, u1))*u1;
u3=T1(3,:)'-(dot(u1, T1(3,:)')/dot(u1, u1))*u1-(dot(u2, T1(3,:)')/dot(u2, u2))*u2;
u4=T1(4,:)'-(dot(u1, T1(4,:)')/dot(u1, u1))*u1-(dot(u2, T1(4,:)')/dot(u2, u2))*u2-(dot(u3, T1(4,:)')/dot(u3, u3))*u3;
u5=T1(5,:)'-(dot(u1, T1(5,:)')/dot(u1, u1))*u1-(dot(u2, T1(5,:)')/dot(u2, u2))*u2-(dot(u3, T1(5,:)')/dot(u3, u3))*u3-(dot(u4, T1(5,:)')/dot(u4, u4))*u4;

v1=dot(v,u2)/dot(u2,u2)*u2+dot(v,u3)/dot(u3,u3)*u3+dot(v,u4)/dot(u4,u4)*u4+dot(v,u5)/dot(u5,u5)*u5;
v2=v-u1-v1;
%%%%START HERE MARTIN
% %Normalized projection of given vector v1 component into M1 space
p1=dot(v1,u2)/dot(u2, u2)*u2+dot(v1,u3)/dot(u3, u3)*u3+dot(v1,u4)/dot(u4, u4)*u4+dot(v1,u5)/dot(u5, u5)*u5;
n1=(1/norm(p1))*p1;
% %Normalized projection of T1 rows into M1 space
f1=dot(T1(1,:)',u2)/dot(u2, u2)*u2+dot(T1(1,:)',u3)/dot(u3, u3)*u3+dot(T1(1,:)',u4)/dot(u4, u4)*u4+dot(T1(1,:)',u5)/dot(u5, u5)*u5;
g1=(1/norm(f1))*f1;
f2=dot(T1(2,:)',u2)/dot(u2, u2)*u2+dot(T1(2,:)',u3)/dot(u3, u3)*u3+dot(T1(2,:)',u4)/dot(u4, u4)*u4+dot(T1(2,:)',u5)/dot(u5, u5)*u5;
g2=(1/norm(f2))*f2;
f3=dot(T1(3,:)',u2)/dot(u2, u2)*u2+dot(T1(3,:)',u3)/dot(u3, u3)*u3+dot(T1(3,:)',u4)/dot(u4, u4)*u4+dot(T1(3,:)',u5)/dot(u5, u5)*u5;
g3=(1/norm(f3))*f3;
f4=dot(T1(4,:)',u2)/dot(u2, u2)*u2+dot(T1(4,:)',u3)/dot(u3, u3)*u3+dot(T1(4,:)',u4)/dot(u4, u4)*u4+dot(T1(4,:)',u5)/dot(u5, u5)*u5;
g4=(1/norm(f4))*f4;
f5=dot(T1(5,:)',u2)/dot(u2, u2)*u2+dot(T1(5,:)',u3)/dot(u3, u3)*u3+dot(T1(5,:)',u4)/dot(u4, u4)*u4+dot(T1(5,:)',u5)/dot(u5, u5)*u5;
g5=(1/norm(f5))*f5;
% %Dot products of normalized projections for FIRST EFFECTS SPACE
d1=dot(n1,g1); d2=dot(n1,g2); d3=dot(n1,g3); d4=dot(n1,g4); d5=dot(n1,g5)
% 
% %Calculations for SECOND EFFECTS SPACE
 %u5=T2(5,:)'-(dot(u1, T2(5,:)')/dot(u1, u1))*u1-(dot(u2, T2(5,:)')/dot(u2, u2))*u2-(dot(u3, T2(5,:)')/dot(u3, u3))*u3-(dot(u4, T2(5,:)')/dot(u4, u4))*u4;
 u6=T2(6,:)'-(dot(u1, T2(6,:)')/dot(u1, u1))*u1-(dot(u2, T2(6,:)')/dot(u2, u2))*u2-(dot(u3, T2(6,:)')/dot(u3, u3))*u3-(dot(u4, T2(6,:)')/dot(u4, u4))*u4-(dot(u5, T2(6,:)')/dot(u5, u5))*u5;
 u7=T2(7,:)'-(dot(u1, T2(7,:)')/dot(u1, u1))*u1-(dot(u2, T2(7,:)')/dot(u2, u2))*u2-(dot(u3, T2(7,:)')/dot(u3, u3))*u3-(dot(u4, T2(7,:)')/dot(u4, u4))*u4-(dot(u5, T2(7,:)')/dot(u5, u5))*u5-(dot(u6, T2(7,:)')/dot(u6, u6))*u6; 
 u8=T2(8,:)'-(dot(u1, T2(8,:)')/dot(u1, u1))*u1-(dot(u2, T2(8,:)')/dot(u2, u2))*u2-(dot(u3, T2(8,:)')/dot(u3, u3))*u3-(dot(u4, T2(8,:)')/dot(u4, u4))*u4-(dot(u5, T2(8,:)')/dot(u5, u5))*u5-(dot(u6, T2(8,:)')/dot(u6, u6))*u6-(dot(u7, T2(8,:)')/dot(u7, u7))*u7;
 u9=T2(9,:)'-(dot(u1, T2(9,:)')/dot(u1, u1))*u1-(dot(u2, T2(9,:)')/dot(u2, u2))*u2-(dot(u3, T2(9,:)')/dot(u3, u3))*u3-(dot(u4, T2(9,:)')/dot(u4, u4))*u4-(dot(u5, T2(9,:)')/dot(u5, u5))*u5-(dot(u6, T2(9,:)')/dot(u6, u6))*u6-(dot(u7, T2(9,:)')/dot(u7, u7))*u7-(dot(u8, T2(9,:)')/dot(u8, u8))*u8;
 u10=T2(10,:)'-(dot(u1, T2(10,:)')/dot(u1, u1))*u1-(dot(u2, T2(10,:)')/dot(u2, u2))*u2-(dot(u3, T2(10,:)')/dot(u3, u3))*u3-(dot(u4, T2(10,:)')/dot(u4, u4))*u4-(dot(u5, T2(10,:)')/dot(u5, u5))*u5-(dot(u6, T2(10,:)')/dot(u6, u6))*u6-(dot(u7, T2(10,:)')/dot(u7, u7))*u7-(dot(u8, T2(10,:)')/dot(u8, u8))*u8-(dot(u9, T2(10,:)')/dot(u9, u9))*u9;
 %u10=T2(5,:)'-(dot(u1, T2(5,:)')/dot(u1, u1))*u1-(dot(u2, T2(5,:)')/dot(u2, u2))*u2-(dot(u3, T2(5,:)')/dot(u3, u3))*u3-(dot(u4, T2(5,:)')/dot(u4, u4))*u4-(dot(u5, T2(5,:)')/dot(u5, u5))*u5-(dot(u6, T2(5,:)')/dot(u6, u6))*u6-(dot(u7, T2(5,:)')/dot(u7, u7))*u7-(dot(u8, T2(5,:)')/dot(u8, u8))*u8-(dot(u9, T2(5,:)')/dot(u9, u9))*u9;
 p2=dot(v2,u6)/dot(u6, u6)*u6+dot(v2,u7)/dot(u7, u7)*u7+dot(v2,u8)/dot(u8, u8)*u8+dot(v2,u9)/dot(u9, u9)*u9+dot(v2,u10)/dot(u10, u10)*u10;
 n2=(1/norm(p2))*p2;
% %Normalized projection of T2 rows into M2 space
 j1=dot(T2(1,:)',u6)/dot(u6, u6)*u6+dot(T2(1,:)',u7)/dot(u7, u7)*u7+dot(T2(1,:)',u8)/dot(u8, u8)*u8+dot(T2(1,:)',u9)/dot(u9, u9)*u9+dot(T2(1,:)',u10)/dot(u10, u10)*u10;
 k1=(1/norm(j1))*j1;
 j2=dot(T2(2,:)',u6)/dot(u6, u6)*u6+dot(T2(2,:)',u7)/dot(u7, u7)*u7+dot(T2(2,:)',u8)/dot(u8, u8)*u8+dot(T2(2,:)',u9)/dot(u9, u9)*u9+dot(T2(2,:)',u10)/dot(u10, u10)*u10;
 k2=(1/norm(j2))*j2;
 j3=dot(T2(3,:)',u6)/dot(u6, u6)*u6+dot(T2(3,:)',u7)/dot(u7, u7)*u7+dot(T2(3,:)',u8)/dot(u8, u8)*u8+dot(T2(3,:)',u9)/dot(u9, u9)*u9+dot(T2(3,:)',u10)/dot(u10, u10)*u10;
 k3=(1/norm(j3))*j3;
 j4=dot(T2(4,:)',u6)/dot(u6, u6)*u6+dot(T2(4,:)',u7)/dot(u7, u7)*u7+dot(T2(4,:)',u8)/dot(u8, u8)*u8+dot(T2(4,:)',u9)/dot(u9, u9)*u9+dot(T2(4,:)',u10)/dot(u10, u10)*u10;
 k4=(1/norm(j4))*j4;
 j5=dot(T2(5,:)',u5)/dot(u5, u5)*u5+dot(T2(5,:)',u6)/dot(u6, u6)*u6+dot(T2(5,:)',u8)/dot(u8, u8)*u8+dot(T2(5,:)',u9)/dot(u9, u9)*u9+dot(T2(5,:)',u10)/dot(u10, u10)*u10;
 k5=(1/norm(j5))*j5;
 j6=dot(T2(6,:)',u5)/dot(u5, u5)*u5+dot(T2(6,:)',u6)/dot(u6, u6)*u6+dot(T2(6,:)',u8)/dot(u8, u8)*u8+dot(T2(6,:)',u9)/dot(u9, u9)*u9+dot(T2(6,:)',u10)/dot(u10, u10)*u10;
 k6=(1/norm(j6))*j6;
 j7=dot(T2(7,:)',u5)/dot(u5, u5)*u5+dot(T2(7,:)',u6)/dot(u6, u6)*u6+dot(T2(7,:)',u8)/dot(u8, u8)*u8+dot(T2(7,:)',u9)/dot(u9, u9)*u9+dot(T2(7,:)',u10)/dot(u10, u10)*u10;
 k7=(1/norm(j7))*j7;
 j8=dot(T2(8,:)',u5)/dot(u5, u5)*u5+dot(T2(8,:)',u6)/dot(u6, u6)*u6+dot(T2(8,:)',u8)/dot(u8, u8)*u8+dot(T2(8,:)',u9)/dot(u9, u9)*u9+dot(T2(8,:)',u10)/dot(u10, u10)*u10;
 k8=(1/norm(j8))*j8;
 j9=dot(T2(9,:)',u5)/dot(u5, u5)*u5+dot(T2(9,:)',u6)/dot(u6, u6)*u6+dot(T2(9,:)',u8)/dot(u8, u8)*u8+dot(T2(9,:)',u9)/dot(u9, u9)*u9+dot(T2(9,:)',u10)/dot(u10, u10)*u10;
 k9=(1/norm(j9))*j9;
 j10=dot(T2(10,:)',u5)/dot(u5, u5)*u5+dot(T2(10,:)',u6)/dot(u6, u6)*u6+dot(T2(10,:)',u8)/dot(u8, u8)*u8+dot(T2(10,:)',u9)/dot(u9, u9)*u9+dot(T2(10,:)',u10)/dot(u10, u10)*u10;
 k10=(1/norm(j10))*j10;
%%Dot products of normalized projections for SECOND EFFECTS SPACE
 e1=dot(n2,k1); e2=dot(n2,k2); e3=dot(n2,k3); e4=dot(n2,k4); e5=dot(n2,k5); e6=dot(n2,k6); e7=dot(n2,k7); e8=dot(n2,k8); e9=dot(n2,k9); e10=dot(n2,k10);
% 

%Displayed Information
fprintf('\n The zero effect vector is \n') 
disp (u1);
fprintf('\n The first effect vector is \n')
disp (v1);
fprintf('\n The second effect vector is \n')
disp (v2);
    fprintf('\n The alignment values of the FIRST effects space are \n')
disp('A B C D E')
disp ([d1 d2 d3 d4 d5]);
    fprintf('\n The alignment values of the SECOND effects space are \n')
disp('AB AC AD AE BC BD BE CD CE DE')
disp ([e1 e2 e3 e4 e5 e6 e7 e8 e9 e10]);