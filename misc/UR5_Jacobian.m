clear all
close all

% UR kinematics parameters
syms t1 t2 t3 t4 t5 t6 real; % theta1 - theta6
syms a2 a3 d1 d4 d5 d6 real;  % DH parameters
syms s1 c1 s2 c2 real;  % sin cos

s1=sin(t1); c1=cos(t1);

s2=sin(t2); c2=cos(t2);

s3=sin(t3); c3=cos(t3);

s4=sin(t4); c4=cos(t4);

s5=sin(t5); c5=cos(t5);

s6=sin(t6); c6=cos(t6);

%UR5 parameters, m
a2 = -0.425; a3 = -0.392;
d1 = 0.0895; d4 = 0.10915; d5 = 0.09465; d6 = 0.0823; 

% UR transformation matrices 
H1  = [c1 0  s1 0;...
       s1 0 -c1 0;...
       0  1  0  d1;...
       0  0  0  1];
T01 = H1;   
   
H2 = [c2 -s2 0 a2*c2;...
      s2 c2 0 a2*s2;...
      0  0  1 0;...
      0  0  0 1];
T02 = H1*H2;  
  
H3 = [c3 -s3 0 a3*c3;...
      s3  c3 0 a3*s3;...
      0  0  1  0;...
      0  0  0  1];
T03 = H1*H2*H3;   
  
H4 = [c4 0 s4 0;...
      s4 0 -c4 0;...
      0  1  0  d4;...
      0  0  0  1];
T04 = H1*H2*H3*H4;  
  
H5 = [c5 0 -s5  0;...
      s5 0  c5  0;...
      0 -1  0  d5;...
      0  0  0   1];
T05 = H1*H2*H3*H4*H5;  
  
H6 = [c6 -s6 0 0;...
      s6 c6  0 0;...
      0   0  1 d6;...
      0   0  0 1];
T06 = H1*H2*H3*H4*H5*H6;

% rotation matrixes  
R01 = H1(1:3,1:3); 
R12 = H2(1:3,1:3);   
R23 = H3(1:3,1:3);  
R34 = H4(1:3,1:3);   
R45 = H5(1:3,1:3); 
R56 = H6(1:3,1:3);

% position vectors p
p6 = T06(1:3,4);

% Geometric Jacobain 
% directions of the joint axes z(i-1)
k = [0 0 1]';
z0 = k;
z1 = R01*k;
z2 = R01*R12*k;
z3 = R01*R12*R23*k;
z4 = R01*R12*R23*R34*k;
z5 = R01*R12*R23*R34*R45*k;
  
J = simplify([diff(p6,t1) diff(p6,t2) diff(p6,t3) diff(p6,t4) diff(p6,t5) diff(p6,t6);...
      z0 z1 z2 z3 z4 z5])

% Inverse Jacobian
%inv(J)

% Analytic Jacobian
syms phi theta psi % Euler angles
sphi = sin(phi);   cphi = cos(phi);
st = sin(theta); ct = cos(theta);
spsi = sin(psi);   cpsi = cos(psi);

Ba = [cpsi*st -spsi 0;...
      spsi*st  cpsi 0;...
      cphi     0    1];

I = eye(3);
O = zeros(3);
Ja = [I  O;...
      O inv(Ba)]*J
 
  
invJa = inv(Ja)
  
