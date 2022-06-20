function Spec4(v)
% Spec4(v)
%   
%   Input a column vector from R^6 v=[a1;a2;a3;a4;a5;a6]
%   Outputs three vectors v0, v1, v2 so that 
%   v=v0+v1+v2
%   with each vi an effects vector
%
%   Example: Find 0, 1, and 2 effects vector for [12;3;8;7;16;2]
%
%   v=[12;3;8;7;16;2]
%   Spec4(v) gives v0=[8;8;8;8;8;8]; 
%   v1=[5;-6.5;0.5;-0.5;6.5;-5]; 
%   v2=[-1;-1.5;-0.5;-0.5;1.5;-1]
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
T1=[1 1 1 0 0 0;1 0 0 1 1 0; 0 1 0 1 0 1; 0 0 1 0 1 1];
T2=eye(6);
z0=(1/6)*ones(6,1);
%Calculation of Effects vectors using projections
u1=dot(v,z0)/dot(z0,z0)*z0;
u2=T1(2,:)'-(dot(u1, T1(2,:)')/dot(u1, u1))*u1;
u3=T1(3,:)'-(dot(u1, T1(3,:)')/dot(u1, u1))*u1-(dot(u2, T1(3,:)')/dot(u2, u2))*u2;
u4=T1(4,:)'-(dot(u1, T1(4,:)')/dot(u1, u1))*u1-(dot(u2, T1(4,:)')/dot(u2, u2))*u2-(dot(u3, T1(4,:)')/dot(u3, u3))*u3;
v1=dot(v,u2)/dot(u2,u2)*u2+dot(v,u3)/dot(u3,u3)*u3+dot(v,u4)/dot(u4,u4)*u4;
v2=v-u1-v1;
%Normalized projection of given vector v1 component into M1 space
p1=dot(v1,u2)/dot(u2, u2)*u2+dot(v1,u3)/dot(u3, u3)*u3+dot(v1,u4)/dot(u4, u4)*u4;
n1=(1/norm(p1))*p1;
%Normalized projection of T1 rows into M1 space
f1=dot(T1(1,:)',u2)/dot(u2, u2)*u2+dot(T1(1,:)',u3)/dot(u3, u3)*u3+dot(T1(1,:)',u4)/dot(u4, u4)*u4;
g1=(1/norm(f1))*f1;
f2=dot(T1(2,:)',u2)/dot(u2, u2)*u2+dot(T1(2,:)',u3)/dot(u3, u3)*u3+dot(T1(2,:)',u4)/dot(u4, u4)*u4;
g2=(1/norm(f2))*f2;
f3=dot(T1(3,:)',u2)/dot(u2, u2)*u2+dot(T1(3,:)',u3)/dot(u3, u3)*u3+dot(T1(3,:)',u4)/dot(u4, u4)*u4;
g3=(1/norm(f3))*f3;
f4=dot(T1(4,:)',u2)/dot(u2, u2)*u2+dot(T1(4,:)',u3)/dot(u3, u3)*u3+dot(T1(4,:)',u4)/dot(u4, u4)*u4;
g4=(1/norm(f4))*f4;
%Dot products of normalized projections for FIRST EFFECTS SPACE
d1=dot(n1,g1); d2=dot(n1,g2); d3=dot(n1,g3); d4=dot(n1,g4);

%Calculations for SECOND EFFECTS SPACE
u5=T2(5,:)'-(dot(u1, T2(5,:)')/dot(u1, u1))*u1-(dot(u2, T2(5,:)')/dot(u2, u2))*u2-(dot(u3, T2(5,:)')/dot(u3, u3))*u3-(dot(u4, T2(5,:)')/dot(u4, u4))*u4;
u6=T2(6,:)'-(dot(u1, T2(6,:)')/dot(u1, u1))*u1-(dot(u2, T2(6,:)')/dot(u2, u2))*u2-(dot(u3, T2(6,:)')/dot(u3, u3))*u3-(dot(u4, T2(6,:)')/dot(u4, u4))*u4-(dot(u5, T2(6,:)')/dot(u5, u5))*u5;
%Normalized projection of given vector v2 component into M2 space
p2=dot(v2,u5)/dot(u5, u5)*u5+dot(v2,u6)/dot(u6, u6)*u6;
n2=(1/norm(p2))*p2;
%Normalized projection of T2 rows into M2 space
j1=dot(T2(1,:)',u5)/dot(u5, u5)*u5+dot(T2(1,:)',u6)/dot(u6, u6)*u6;
k1=(1/norm(j1))*j1;
j2=dot(T2(2,:)',u5)/dot(u5, u5)*u5+dot(T2(2,:)',u6)/dot(u6, u6)*u6;
k2=(1/norm(j2))*j2;
j3=dot(T2(3,:)',u5)/dot(u5, u5)*u5+dot(T2(3,:)',u6)/dot(u6, u6)*u6;
k3=(1/norm(j3))*j3;
j4=dot(T2(4,:)',u5)/dot(u5, u5)*u5+dot(T2(4,:)',u6)/dot(u6, u6)*u6;
k4=(1/norm(j4))*j4;
j5=dot(T2(5,:)',u5)/dot(u5, u5)*u5+dot(T2(5,:)',u6)/dot(u6, u6)*u6;
k5=(1/norm(j5))*j5;
j6=dot(T2(6,:)',u5)/dot(u5, u5)*u5+dot(T2(6,:)',u6)/dot(u6, u6)*u6;
k6=(1/norm(j6))*j6;
%Dot products of normalized projections for SECOND EFFECTS SPACE
e1=dot(n2,k1); e2=dot(n2,k2); e3=dot(n2,k3); e4=dot(n2,k4); e5=dot(n2,k5); e6=dot(n2,k6);


%Displayed Information
fprintf('\n The zero effect vector is \n') 
disp (u1);
fprintf('\n The first effect vector is \n')
disp (v1);
fprintf('\n The second effect vector is \n')
disp (v2);
fprintf('\n The alignment values of the FIRST effects space are \n')
disp('G H R S')
disp ([d1 d2 d3 d4]);
fprintf('\n The alignment values of the SECOND effects space are \n')
disp('GH GR GS HR HS RS')
disp ([e1 e2 e3 e4 e5 e6]);