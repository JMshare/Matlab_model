%% To convert the A_body in \dot{v}=A_body.v (v=[u,v,w]) to A_aero                  as A_aero = A_body.inv(R)
%%         or the A_wind in \dot{wind}=A_wind.wind (wind=[V,Alpha,Beta]) to A_body  as A_body = A_wind.R

%% Get R_inv = v'V and v'Alpha and v'Beta

V = 5.199999607713583;
B = 0;
A = deg2rad(45.891366219060579);
du = 0.001;
Apl = A + du;
Amn = A - du;
[upl,vpl,wpl] = Extract_body_states(V, rad2deg(Apl), rad2deg(B));
[umn,vmn,wmn] = Extract_body_states(V, rad2deg(Amn), rad2deg(B));
du_Alp = (upl - umn)/(2*du);
dv_Alp = (vpl - vmn)/(2*du);
dw_Alp = (wpl - wmn)/(2*du);
dv_Alp = [du_Alp; dv_Alp; dw_Alp];

B = 0;
A = deg2rad(45.891366219060579);
du = 0.001;
Bpl = B + du;
Bmn = B - du;
[upl,vpl,wpl] = Extract_body_states(V, rad2deg(A), rad2deg(Bpl));
[umn,vmn,wmn] = Extract_body_states(V, rad2deg(A), rad2deg(Bmn));
du_Bt = (upl - umn)/(2*du);
dv_Bt = (vpl - vmn)/(2*du);
dw_Bt = (wpl - wmn)/(2*du);
dv_Bt = [du_Bt; dv_Bt; dw_Bt];

B = 0;
A = deg2rad(45.891366219060579);
du = 0.001;
Vpl = V + du;
Vmn = V - du;
[upl,vpl,wpl] = Extract_body_states(Vpl, rad2deg(A), rad2deg(B));
[umn,vmn,wmn] = Extract_body_states(Vmn, rad2deg(A), rad2deg(B));
du_V = (upl - umn)/(2*du);
dv_V = (vpl - vmn)/(2*du);
dw_V = (wpl - wmn)/(2*du);
dv_V = [du_V; dv_V; dw_V];

R_inv = [dv_V, dv_Alp, dv_Bt];


%% Get R = V'v and Alpha'v and Beta'v
[u,v,w] = Extract_body_states(V, rad2deg(A), rad2deg(B));
v = [u, v, w];
vp = [du, 0, 0];
dAlpha_u = (Alpha(v+vp) - Alpha(v-vp))/(2*du) % -u^-2*w / ((w/u)^2 + 1)
-v(3) / ((v(1)^2)*((v(3)/v(1))^2 + 1))

vp = [0, 0, du];
dAlpha_w = (Alpha(v+vp) - Alpha(v-vp))/(2*du) % u^-1 / ((w/u)^2 + 1)
1 / (v(1)*((v(3)/v(1))^2 + 1))

vp = [0, du, 0];
dAlpha_v = (Alpha(v+vp) - Alpha(v-vp))/(2*du) % u^-1 / ((w/u)^2 + 1)
0


vp = [du, 0, 0];
dBeta_u = (Beta(v+vp) - Beta(v-vp))/(2*du)
0 % iff vN == 0

vp = [0, du, 0];
dBeta_v = (Beta(v+vp) - Beta(v-vp))/(2*du)
1/norm(v) % iff vN == 0

vp = [0, 0, du];
dBeta_w = (Beta(v+vp) - Beta(v-vp))/(2*du)
0 % iff vN == 0


vp = [du, 0, 0];
dV_u = (norm(v+vp) - norm(v-vp))/(2*du)
v(1)/norm(v)

vp = [0, du, 0];
dV_v = (norm(v+vp) - norm(v-vp))/(2*du)
v(2)/norm(v) 

vp = [0, 0, du];
dV_w = (norm(v+vp) - norm(v-vp))/(2*du)
v(3)/norm(v) 

R = [dV_u, dV_v, dV_w;  dAlpha_u, dAlpha_v, dAlpha_w; dBeta_u, dBeta_v, dBeta_w];


%% Compare the R a R_inv
[R_inv, inv(R)]


%% Load A from Run_A_eigs at the V specified here and compare
[A(7:9, 4:6)*R_inv, fAero(7:9,:)]
[fAero(7:9,:)*R, A(7:9,4:6)]




function [alpha] = Alpha(v)
    alpha = atan(v(3)/v(1));
end


function [beta] = Beta(v)
    beta = asin(v(2)/norm(v));
end

