function Spec4(v)
% Spec4(v)
%   
%   Input a column vector from R^6 v=[a1;a2;a3;a4;a5;a6]
%   Outputs three vectors v0, v1, v2 so that 
%   v=v0+v1+v2
%   with each vi an effects vector
%
%   Example: Find 0, 1, and 2 effects vector for [12;3;8;7;16;2]
%   v=[12;3;8;7;16;2]
%   Spec4(v) gives v0=[8;8;8;8;8;8] ; v1=[5;-6.5;0.5;-0.5;6.5;-5]; v2=[-1;-1.5;-0.5;-0.5;1.5;-1]

T1=[1 1 1 0 0 0;1 0 0 1 1 0; 0 1 0 1 0 1; 0 0 1 0 1 1];
T2=eye(6);
z0=(1/6)*ones(6,1);
u1=dot(v,z0)/dot(z0,z0)*z0;
u2=T1(2,:)'-(dot(u1, T1(2,:)')/dot(u1, u1))*u1;
u3=T1(3,:)'-(dot(u1, T1(3,:)')/dot(u1, u1))*u1-(dot(u2, T1(3,:)')/dot(u2, u2))*u2;
u4=T1(4,:)'-(dot(u1, T1(4,:)')/dot(u1, u1))*u1-(dot(u2, T1(4,:)')/dot(u2, u2))*u2-(dot(u3, T1(4,:)')/dot(u3, u3))*u3;
v1=dot(v,u2)/dot(u2,u2)*u2+dot(v,u3)/dot(u3,u3)*u3+dot(v,u4)/dot(u4,u4)*u4;
v2=v-u1-v1;
fprintf('\n The zero effect vector is \n') 
disp (u1);
fprintf('\n The first effect vector is \n')
disp (v1);
fprintf('\n The second effect vector is \n')
disp (v2);