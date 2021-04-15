%%% Comparing my calculations with Etkin p.166 %%%

cmq = -23.92;
clp = -0.3295;
cmde = -1.444;
clda = -1.368*1e-2;
rho = 0.3045;
S = 511;
V = 235.9;
b = 59.64;
c = 8.324;
mq = 0.25*rho*V*c^2*S*cmq;
lp = 0.25*rho*V*b^2*S*clp;
mde = 0.5*rho*(V^2)*S*c*cmde;
lda = 0.5*rho*(V^2)*S*b*clda;

mq_star = -1.521*1e7;
lp_star = -1.076*1e7;
mde_star = -3.839*1e7*1.355818; % to Nm from ftlb
lda_star = 0;


[mq, mq_star]
[lp, lp_star]
[mde, mde_star]



%% A matrix
Ix = 0.247*1e8;
Iy = 0.449*1e8;
Iz = 0.673*1e8;
Izx = -0.212*1e7;
I = [Ix, 0, Izx; 0, Iy, 0; Izx, 0, Iz];
m = (2.83176*1e6)/9.81;

mwdot_star = -1.702*1e4;
zq_star = -4.524*1e5;
zwdot_star = 1.909*1e3;

Aqdot = (mq_star + mwdot_star*(zq_star+m*V)/(m-zwdot_star))/Iy;
Aqdot_star = -0.4285;
[Aqdot, Aqdot_star]




