function qddot = fun_qddotsspsinglecontact(x,u)

global A B Nx Nu pert MI L m  nx ny tx ty g r lam vars misc alp alpval indic kc lamall xdata lamx lamy val af acal fx fy Mmat2 invM phi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1 = m(2);
m2 = m(3);

L1 = L(1);
L2 = L(2);

r1 = r(2);
r2 = r(3);

MI1 = MI(2);
MI2 = MI(3);

%%%%%%%%%%%%%%%%%%%%%%
tht1 = x(3);
tht2 = x(4);
xh   = x(1);
yh   = x(2);
omg1 = x(7);
omg2 = x(8);
vhx = x(5);
vhy = x(6);


F1 = u(3);
F2 = u(4);
T1  =  u(1);
T2  =  u(2);

c = -1.6;
%T1  = u(1);
%{
T2  = c*omg2*r2 +u(1);
T3  = c*omg3*r3 +u(2);
T4  = c*omg4*r4 +u(3);
T5  = c*omg5*r5 +u(4);
T6  = c*omg6*r6 +u(5);
T7  = c*omg7*r7 +u(6);
%}
lam1 = lamy;
lam2 =lamx;
%{
T2  = -c*omg2*r2 +u(1);
T3  = -c*omg3*r3 +u(2);
T4  = -MI4*alp4  -c*omg4*r4 +u(3);
T5  = -c*omg5*r5 +u(4);
T6  = -c*omg6*r6 +u(5);
T7  = -c*omg7*r7 +u(6);

%}

%T4  =  c*omg4 ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mmat=[m1+m2,0,L1.*m2.*sin(tht1)+m1.*r1.*sin(tht1),m2.*r2.*sin(tht2);0, ...
  m1+m2,(-1).*L1.*m2.*cos(tht1)+(-1).*m1.*r1.*cos(tht1),(-1).*m2.* ...
  r2.*cos(tht2);L1.*m2.*sin(tht1)+m1.*r1.*sin(tht1),(-1).*L1.*m2.* ...
  cos(tht1)+(-1).*m1.*r1.*cos(tht1),L1.^2.*m2+MI1+m1.*r1.^2,L1.*m2.* ...
  r2.*cos(tht1+(-1).*tht2);m2.*r2.*sin(tht2),(-1).*m2.*r2.*cos(tht2) ...
  ,L1.*m2.*r2.*cos(tht1+(-1).*tht2),MI2+m2.*r2.^2];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mmat'
Mmat - Mmat'
issymmetric(Mmat)
eig(Mmat)


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% phi is for finding qddot with ext lambda calculation 
phi=[F1+(-1).*L1.*m2.*omg1.^2.*cos(tht1)+(-1).*m1.*omg1.^2.*r1.*cos( ...
  tht1)+(-1).*m2.*omg2.^2.*r2.*cos(tht2),F2+(-1).*g.*m1+(-1).*g.*m2+ ...
  (-1).*L1.*m2.*omg1.^2.*sin(tht1)+(-1).*m1.*omg1.^2.*r1.*sin(tht1)+ ...
  (-1).*m2.*omg2.^2.*r2.*sin(tht2),T1+(-1).*T2+g.*L1.*m2.*cos(tht1)+ ...
  g.*m1.*r1.*cos(tht1)+(-1).*L1.*m2.*omg2.^2.*r2.*sin(tht1+(-1).* ...
  tht2),T2+g.*m2.*r2.*cos(tht2)+L1.*m2.*omg1.^2.*r2.*sin(tht1+(-1).* ...
  tht2)];



qddot_invdy = invM*phi'
af = qddot_invdy ;
qddot = qddot_invdy;
 

end






 






