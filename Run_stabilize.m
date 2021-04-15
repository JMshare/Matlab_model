%% Run LQR with Integral and Filter (control time delay) and Kalman filter 
% K stabilization from LQR for states = [z,y,z, u,v,w, p,q,r, phi,theta,psi]
% Integral hold for [x,y,z, phi,theta,psi]
% 1st order lag effect for control actuation acting as low pass filter on
% control
% Kalman filter for measurement and process noise
% can use nonlinear ode_fun if filter on
% last non-experimental version after 16

clear;
close all;
clc;
format;

global y_pert_e

vars_name = 'vars_fits';

%% Setup
control = 1;

bin_int = []; % stabilization only
% bin_int = [3, 12]; % altitude cntrl
% bin_int = [1:3, 12]; % pos cntrl
% bin_int = [10,11,12]; % att control, possible only with filter control on, use when position sensors not reliable and set binr = [123456]
% bin_int = [10,11];

filter = 0; % design control such that it includes for the lags/filters
tf = 100./[1; 1; 1; 1]; % time constants for control lag filter
delay = 0; % simulate delay in control

binr = [1,2,3,4,5,6]; % remove from feedback control stabilization
% binr = [];

estim = 0;
nonlin = 0;

design = 0;
% w = [1,1,1, 1,1,1,10, 1, 1, 1, 100]; % [binx, binu, bine,  nux, nuv, nuomg, nueps, mp, mq, mr,  mu] % weights in LQR obj fn J
% w = [1,1,60, 1,1,1,46.1498, 1, 1.25, 77.49, 111]; % [binx, binu, bine,  nux, nuv, nuomg, nueps, mp, mq, mr,  mu] % weights in LQR obj fn J
w = [1,1,60, 1,1,1,50, 1, 1, 1, 100];
MC = 1000;

% Fs = [2:2:20, 35, 50, 100, 250, 300]; % discrete freqns for which to plot eigs
Fs = [10, 5, 10, 20]; % it will assume discrete-time system if the first element is > 0

%% Inputs
Npts = 15000;
T = 5;
Tstep = 1;
y_pert_v = 1*[0; 0; 0; 0;0;0; 0;0;0; 0;0;0];
y_pert_e = 1*[0; 0; 0; 0;0;0; 0.0;0.0;0.0; deg2rad(45);deg2rad(0);deg2rad(0)];
y_pert = 1*y_pert_v + 1*y_pert_e;

y_step = 0*[0;0;0; 0;0;0; 0;0;0; 0.0;0.0;0.0];
u_step = 0*[1;1;1;1];

proc_bias = 0*[0;0;0; 1;1;1; 0;0;0; 0;0;0];

%% Checks
bin_int = bin_int(~ismember(bin_int, binr)); % remove integral if removed from stabilization

if (filter == 1) && (delay == 0) % if control with filter, then simulate it
	delay = 1;
	fprintf('Filter design without simulating it is not supported!\n Running with simulation of the lags/filters\n');
end

if filter == 0 && nonlin == 1
    nonlin = 0;
    fprintf('Nonlinear case when filter is off is not supported!\n Running for nonlin=0\n');
end

%% Loads
load(sprintf('%s_syst.mat', vars_name), 'A', 'Argrbp', 'B', 'Brgrbp', 'Data');
yN = Data.yN;
uN = Data.uN;

A = Argrbp;
A(4:6, 1:6) = 0; % I dont have a control over x,v so they will make it unstable whatever i do, plus this doesnt simulate the real instability I am getting cos my v is essentially zero, it wiggles on a spot.

A = 0*A;
% A(1:3, 10:12) = [0, -9.81, 0; 9.81, 0, 0; 0, 0, 0]; % add g

A = [[zeros(3), eye(3), zeros(3,6)]; A; [zeros(3,6), eye(3), zeros(3)]];
%B = Brgrbp;
B = sign(B).*mean(abs(B),2);
B = B;
B(1:3,:) = 0; % wont simulate the u,v,w but helps if i want to look at the final eigenvalues
B = [zeros(3,size(B,2)); B; zeros(3,size(B,2))];
y0 = 0*yN(1,:)';
u0 = mean(uN, 1)';
u0 = (0*u0+1).*mean(u0);
C = eye(length(y0)); % not necessarily if eps0~=0, check the g(y) in data process fn
D = zeros(length(y0), size(B,2));

%% Setpoint step
yd = @(t) y0 + y_step.*(t>Tstep);
ud = @(t) u0 + u_step.*(t>Tstep);


%% Kalman filter dynamics
fprintf('Observability rank %d out of %d on 4th power \n', rank([C; C*A; C*A^2; C*A^3; C*A^4]), length(y0));

G = diag([2;2;2; 2;2;2; 0.05;0.05;0.05; 0.05;0.05;0.05]); % measurement noise
F = diag([0;0;0; 0.1;0.1;0.1; 0.03;0.03;0.03; 0;0;0]); % process noise (nothing will teleport us from one x to other, or change speed in one step => only accelerations noise)
R = G*G';
Q = F*F';
[P, ~, ~] = care(A', C', Q, R);
Kk = P*C'*inv(R);
% [~, Kk, P] = kalman(ss(A,[B eye(12)],C,[D zeros(12)]), Q, R, 0*y0*y0');



%% Controllability
fprintf('Controllability rank %d out of %d on 4th power \n', rank([B, A*B, A^2*B, A^3*B, A^4*B]), length(y0));


%% Integral and Filter control
bini = [1,2,3,10,11,12];
Ci = C(bini,:);
Di = D(bini,:);
r = 0*Ci*y0;

Tf = diag(tf);



%% Assemble the matrices
Z12 = zeros(length(y0), length(y0));
Z21 = Z12;
Z14 = zeros(length(y0), length(bini));
Z24 = Z14;
Z31 = zeros(length(u0), length(y0));
Z32 = Z31;
Z34 = zeros(length(u0), length(bini));
Z41 = zeros(length(bini), length(y0));
Z44 = zeros(length(bini), length(bini));
Zc1 = zeros(length(y0), length(u0));
Zc2 = Zc1;
Zc4 = zeros(length(bini), length(u0));
Zn12 = zeros(length(y0), size(G,2));
Zn21 = zeros(length(y0), size(F,2));
Zn31 = zeros(length(u0), size(F,2));
Zn32 = zeros(length(u0), size(G,2));
Zn41 = zeros(length(bini), size(F,2));
Zn42 = zeros(length(bini), size(G,2));

z0 = [y0; y0; u0; r];
z_pert = [y_pert; y_pert; 0*u0; 0*r];
zd = @(t) [yd(t); yd(t); ud(t); r*0*t];

Aht = [A,    Z12, B, Z14; 
       Z21,  A,   B, Z24;
       Z31,  Z32, -Tf, Z34; 
       Z41,  Ci,  Di, Z44];
Bht = [Zc1; Zc2; Tf; Zc4];
Nht = [F,    Zn12;
       Zn21, Kk*G;
       Zn31, Zn32;
       Zn41, Zn42]; % noise input
Eht = [0*A,  Z12,   0*B,   Z14; 
       Kk*C, -Kk*C, 0*B,   Z24;
       Z31,  Z32,   0*-Tf, Z34; 
       Z41,  0*Ci,  0*Di,  Z44]; % estimator

H = eye(length(z0));
biny = 1:length(y0);
Hy = H(biny, :); % noisy state
binx = length(y0)+1: 2*length(y0);
Hx = H(binx, :); % estimated state
binu = 2*length(y0)+1: 2*length(y0)+length(u0);
Hu = H(binu, :); % filtered control
bine = 2*length(y0)+length(u0)+1:length(z0);
He = H(bine, :); % integral of the error

binlqr = z0==z0; % creating a basket of states used for feedback
binlqr(biny) = 0; % the first y0 states are just for our integration here to get the noisy y. Not used for feedback.
binlqr(bine) = ismember(bini, bin_int); % integral states that we don't want to use for feedback (need to remove here otherwise insufficient rank)
binlqr(binr + length(y0)) = 0; % removing the unwanted states from feedback

%% not Filter/Delay matrices
Aht2 = [A,    Z12, 0*B, Z14; 
        Z21,  A,   0*B, Z24;
        Z31,  Z32, 0*-Tf, Z34; 
        Z41,  Ci,  Di, Z44];
Bht2 = [B; B; 0*Tf; Zc4];
binlqr2 = binlqr;
binlqr2((2*length(y0)+1) : (2*length(y0)+length(u0))) = 0; % remove from feedback matrix cos otherwise insufficient rank


%% LQR design (find the scaling of terms in the Q R N matrices)
bins.binlqr = binlqr;
bins.biny = biny;
bins.binx = binx;
bins.binu = binu;
bins.bine = bine;
bins2 = bins;
bins2.binlqr = binlqr2;
if design ~= 0
    [w] = LQR_design(Aht2, Bht2, w, bins2); 
end
[Q, R, N] = LQR_weight_matrices(size(Aht,1), size(Bht,2), w, bins);


%% Calculate the feedback control  matrix K
Kht = zeros(length(u0),length(z0));
if control ~= 0
	if filter == 1
		[Klqr] = lqr(Aht(binlqr,binlqr), Bht(binlqr,:), Q(binlqr,binlqr), R, N(binlqr,:)); % u = -K.x
		Kht(:,binlqr) = Klqr;
	else
		[Klqr] = lqr(Aht2(binlqr2,binlqr2), Bht2(binlqr2,:), Q(binlqr2,binlqr2), R, N(binlqr2,:)); % u = -K.x
        Kht(:,binlqr2) = Klqr;
	end
end
% p = [  -7.4544 + 3.3330i
%   -7.4544 - 3.3330i
%   -6.8941 + 3.4702i
%   -6.8941 - 3.4702i
%   -2.0184 + 1.8709i
%   -2.0184 - 1.8709i];
% Kpp = place(Aht2(binlqr2,binlqr2), Bht2(binlqr2,:), flipud(p)); % pole-placement
% Kht(:,binlqr2) = Kpp;


%% Monte-Carlo design checks
Monte_Carlo_design(MC, Kht, Aht2, Bht2, bins2, 0);
if Fs(1) > 0
    Monte_Carlo_design(MC, Kht, Aht2, Bht2, bins2, 1/Fs(1));
end


%% Simulation type setup
if (delay == 0) && (filter == 0)
	Aht = Aht2;
	Bht = Bht2;
end


%% Integrate
rng(1);
t = linspace(0,T,Npts)';
dt = t(2) - t(1);
ode_fun = @(y,u) proc_bias; % The dy = f(y,u). If at critocal point, f(y0,u0)=0. If not, problem
if nonlin == 1
    Aht2 = [0*A,    Z12, 0*B, Z14;
            Z21,  0*A,   0*B, Z24;
            Z31,  Z32, -Tf, Z34;
            Z41,  Ci,  Di, Z44]; % only the filter and int states are propagated linearly
    f = @(t,z) [ode_fun(Hy*z,Hu*z)  ; ode_fun(Hx*z,Hu*z)  ; 0*u0 ;r] + Aht2*(z-zd(t)) + Bht*bound_control(-Kht*(z-zd(t)), ud(t)) + estim*(Eht*(z-zd(t)) + Nht*([2*(rand(size(F,2),1)-0.5)/sqrt(dt); 2*(rand(size(G,2),1)-0.5)]));
else
    f = @(t,z) [ode_fun(yd(t),ud(t)); ode_fun(yd(t),ud(t)); 0*u0 ;r] + Aht *(z-zd(t)) + Bht*bound_control(-Kht*(z-zd(t)), ud(t)) + estim*(Eht*(z-zd(t)) + Nht*([2*(rand(size(F,2),1)-0.5)/sqrt(dt); 2*(rand(size(G,2),1)-0.5)]));
end

if Fs(1) > 0
    fs = Fs(1);
    dt = 1/fs;
    t = [0:dt:T]';
    z = ode1d(Aht, Bht, Kht, z0, t, ud)';
else
    z = ode1(f, t, z0);
end



%% Perturb
tp = [T:dt:2*T]';
if Fs(1) > 0
    zp = ode1d(Aht, Bht, Kht, z(end,:)'+z_pert, tp, ud)';
else
    zp = ode1(f, tp, z(end,:)'+z_pert);
end
t = [t(1:end-1); tp];
z = [z(1:end-1, :); zp];


%% Postprocess
z = z';
t = t';
us = bound_control(-Kht*(z-zd(t)), ud(t)) + ud(t);

y = Hy*z; % noisy state
x = Hx*z; % estimated states
uf = Hu*z; % filtered control
eint = He*z; % integral of the error



%% Plots
M = Aht-Bht*Kht;
bina = length(y0)+1 : length(z0);
% ismember(eig(M(bina,bina)), eig(M)) % => can split it like this cos of block matrix properties
plot_eigs(eig(M(bina,bina)), eig(A), [], 10, 0);
plot_eigs(eig(M(binlqr2,binlqr2)), eig(Aht(binlqr2,binlqr2)), [], 11, 0);
plot_trajectories(y, x, t, us, uf, eint, yd(t), 0);
%Animate_trajectory_FG(t, y, 0);
[tau, zeta, Mp] = perform_parameters(eig(M(bina,bina)));
[tau, zeta, Mp]
-Klqr*[0;0;0;deg2rad(0);deg2rad(0);deg2rad(45)]


%% Output for PX4
if any(any([Kht(:,binx), Kht(:,binu), Kht(:,bine)] ~= Kht(:,13:end)))
    fprintf('Error in dimensions!\n');
else
    PX4_printouts(Kht, u0, binx, binu, bine);
end

%% Discrete eigenvalues
for i = 2:length(Fs)
    F = Fs(i);
    sysc = ss(Aht, Bht, eye(size(Aht)), zeros(size(Bht)));
    sysd = c2d(sysc, 1/F);
    Md = sysd.A - sysd.B*Kht;
    sysc0 = ss(Aht,zeros(size(Aht)), zeros(size(Aht)), zeros(size(Aht)));
    sysd0 = c2d(sysc0, 1/F);
    Ahtd = sysd0.A;
    plot_eigs(eig(Md(bina,bina)), eig(Ahtd(bina,bina)), [], 11+i, 1/F);
    plot_eigs(eig(Md(binlqr2,binlqr2)), eig(Ahtd(binlqr2,binlqr2)), [], 21+i, 1/F);
    ps = get(gcf, 'position');
    ps(4) = ps(4)*0.7;
    ps(3) = ps(3)*0.7;
    set(gcf, 'position', ps);
end


%%-------------------------------------------------------------------------
%% Function definitions
function [Q, R, N] = LQR_weight_matrices(n, m, w, bins)
    Q0 = eye(n);
    R0 = eye(m);
    N0 = zeros(n,m);

    Q = Q0;
    Qx = Q(bins.binx, bins.binx);
    Qx(1:3,1:3) = Qx(1:3,1:3)*w(4); % x
    Qx(4:6,4:6) = Qx(4:6,4:6)*w(5); % v
    Qx(7:9,7:9) = Qx(7:9,7:9)*w(6); % omg
    Qx(10:12,10:12) = Qx(10:12,10:12)*w(7); % eps
    Qx([7,10], [7,10]) = Qx([7,10], [7,10])*w(8); % p,phi
    Qx([8,11], [8,11]) = Qx([8,11], [8,11])*w(9); % q,theta
    Qx([9,12], [9,12]) = Qx([9,12], [9,12])*w(10); % r,psi
    Q(bins.binx, bins.binx) = Qx*w(1);
    Q(bins.binu, bins.binu) = Q(bins.binu, bins.binu)*w(2);
    Q(bins.bine, bins.bine) = Q(bins.bine, bins.bine)*w(3);
    
    R = w(end)*R0;
    N = N0;  
end

function [w] = LQR_design(Aht, Bht, w0, bins)
    K_px4 = [[-0.06,  0.06,  0.06, -0.36,  0.36,  0.36];
             [ 0.06, -0.06,  0.06,  0.36, -0.36,  0.36];
             [ 0.06,  0.06, -0.06,  0.36,  0.36, -0.36];
             [-0.06, -0.06, -0.06, -0.36, -0.36, -0.36]];
    binlqr = bins.binlqr;
    ub = 1e12*ones(size(w0));
    lb = zeros(size(w0));
    Aeq = eye(length(w0));
    bineq = [1,2,3,4,5,6,8];  % [binx, binu, bine,  nux, nuv, nuomg, nueps, mp, mq, mr,  mu]
    w = fmincon(@J, w0, [], [], Aeq(bineq,:), w0(bineq), lb, ub);

    function [o] = J(w)
        [Q, R, N] = LQR_weight_matrices(size(Aht,1), size(Bht,2), w, bins);
        try
            K_lqr = lqr(Aht(binlqr,binlqr), Bht(binlqr,:), Q(binlqr,binlqr), R, N(binlqr,:));
            M = Aht(binlqr,binlqr)-Bht(binlqr,:)*K_lqr;
            M_px4 = Aht(binlqr,binlqr)-Bht(binlqr,:)*K_px4;
            [tau, zeta, Mp] = perform_parameters(eig(M));
            [tau_px4, zeta_px4, Mp_px4] = perform_parameters(eig(M_px4));
            [tau, tau_px4]
            [zeta, zeta_px4]
            [Mp, Mp_px4]
            % o = norm(K_px4 - K_lqr);
            o = norm(zeta - 0.9);
            o = o + norm(tau - [0.15; 0.15; 0.15; 0.15; 0.43; 0.43]);
        catch
            fprintf("LQR design optim infeasible weights w\n");
            o = 1e6;
        end
    end
end

function [tau, zeta, Mp] = perform_parameters(lam)
    try
        tau = log(1-0.632)./real(lam);
    catch
        tau = 0*lam;
    end
    
    zeta = -real(lam)./sqrt(real(lam).^2+imag(lam).^2);
    
    try
    Mp = exp((-pi*zeta)./sqrt(1-zeta.^2));
    catch 
        Mp = 0*lam;
    end
end

function [] = Monte_Carlo_design(MC, Kht, Aht, Bht, bins, Ts)
    binlqr = bins.binlqr;
    K_lqr = Kht(:, binlqr);
    lamtds = zeros(MC,nnz(binlqr));
    qs = zeros(MC, 3+3+3*4);
    for i = 1:MC
        qA =  1*[(2*rand)*2, (2*rand)*2, (2*rand-2)*2]*(i>1);
        qB1 = 0*[randn*5, randn*5, randn*5]*(i>1);
        qB2 = 1*[10*(2*rand(1,4)-1); 8*(2*rand(1,4)-1); 3*(2*rand(1,4)-1)]*(i>1);
        Atd1 = zeros(12); 
        Atd1(7:9, 7:9) = diag(qA);
        Btd1 = zeros(12,4); 
        Btd1(7,:) = qB1(1); Btd1(8,:) = qB1(2); Btd1(9,:) = qB1(3);
        Btd2 = zeros(12,4); 
        Btd2(7,:) = qB2(1,:); Btd2(8,:) = qB2(2,:); Btd2(9,:) = qB2(3,:);
        
        Atd = Aht;
        Btd = Bht;
        
        binx = bins.binx;
        Atd(binx,binx) = Atd(binx,binx) + Atd1;
        Btd(binx,:) = Btd(binx,:) + Btd1 + Btd2;
        
        if Ts > 0
            sysc = ss(Atd(binlqr,binlqr), Btd(binlqr,:), eye(size(Atd(binlqr,binlqr))), zeros(size(Btd(binlqr,:))));
            sysd = c2d(sysc, Ts);
            Mtd = sysd.A - sysd.B*K_lqr;
        else
            Mtd = Atd(binlqr,binlqr)-Btd(binlqr,:)*K_lqr;
        end
        
        lamtds(i,:) = eig(Mtd);
        qs(i,:) = [qA, qB1, qB2(:)'];
    end
    
    if MC > 0
        fign = 50 + 10*(Ts>0);
        plot_eigs(lamtds(1,:)', [], lamtds(:), fign, Ts);
        if Ts > 0
            lamtds = log(lamtds)/Ts;
            plot_eigs(lamtds(1,:)', [], lamtds(:), fign+1, 0);
        end
        [tautds, zetatds, Mptds] = perform_parameters(lamtds);
        figure(fign+2)
        histogram(tautds(:), 50, 'Normalization', 'probability');
        title('M-C simulation results, settling times');
        xlabel('\tau_{632} [s]');
        ylabel('probability');
        ps = get(gcf, 'position');
        ps(4) = ps(4)*0.5;
        ps(3) = ps(3)*0.5;
        set(gcf, 'position', ps);
        figure(fign+3)
        histogram(zetatds(:), 50, 'Normalization', 'probability');
        title('M-C simulation results, damping ratios');
        xlabel('\zeta');
        ylabel('probability');
        ps = get(gcf, 'position');
        ps(4) = ps(4)*0.5;
        ps(3) = ps(3)*0.5;
        set(gcf, 'position', ps);
        figure(fign+4)
        histogram(Mptds(:), 50, 'Normalization', 'probability');
        title('M-C simulation results, peak overshoots');
        xlabel('M_p [%]');
        ylabel('probability');
        ps = get(gcf, 'position');
        ps(4) = ps(4)*0.5;
        ps(3) = ps(3)*0.5;
        set(gcf, 'position', ps);
        
%         path = fullfile('..', '..', 'Report','Imgs2'); 
%         if (Ts == 0) && (nnz(qB2) > 0) && (nnz(qA) == 0)
%             filename = fullfile(path, 'LQR_eigs_Btd');
%             saveas(fign, filename, 'epsc');
%             filename = fullfile(path, 'LQR_taus_Btd');
%             saveas(fign+2, filename, 'epsc');
%             filename = fullfile(path, 'LQR_zeta_Btd');
%             saveas(fign+3, filename, 'epsc');
%             filename = fullfile(path, 'LQR_Mp_Btd');
%             saveas(fign+4, filename, 'epsc');
%         end
%         if (Ts == 0) && (nnz(qB2) == 0) && (nnz(qA) > 0)
%             filename = fullfile(path, 'LQR_eigs_Atd');
%             saveas(fign, filename, 'epsc');
%             filename = fullfile(path, 'LQR_taus_Atd');
%             saveas(fign+2, filename, 'epsc');
%             filename = fullfile(path, 'LQR_zeta_Atd');
%             saveas(fign+3, filename, 'epsc');
%             filename = fullfile(path, 'LQR_Mp_Atd');
%             saveas(fign+4, filename, 'epsc');
%         end
%         if (Ts == 0) && (nnz(qB2) > 0) && (nnz(qA) > 0)
%             filename = fullfile(path, 'LQR_eigs_Atd_Btd');
%             saveas(fign, filename, 'epsc');
%         end
    end
end

function [del_u] = bound_control(del_u, u0)
    u_max = [1; 1; 1; 1];
    u_min = [0; 0; 0; 0];
    u = del_u + u0;
    u = max(min(u_max, u), u_min);
    del_u = u - u0;
end


function [] = plot_trajectories(y, x, t, us, uf, eint, yd, binc)


    ttlstr = 'Estimated states';
    plot_states(t, x, 4, ttlstr, yd, binc)
    
    figure(5);
    clf;
    h = plot(t, us);
    set(h, {'DisplayName'}, {'u_1'; 'u_2'; 'u_3'; 'u_4'});
    legend show;
    legend('Location', 'best');
    grid minor;
    title('Control history');
    xlabel('t [s]');
    
    figure(6)
    clf;
    h = plot(t, uf);
    set(h, {'DisplayName'}, {'u_1'; 'u_2'; 'u_3'; 'u_4'});
    legend show;
    legend('Location', 'best');
    grid minor;
    title('Control lagged');
    xlabel('t [s]');
    
    figure(7)
    plot(t, eint);
    grid minor;
    title('Integral error');
    legend('x', 'y', 'z', '\phi', '\theta', '\psi');
    xlabel('t [s]');
    
    
    ttlstr = 'Real states';
    plot_states(t, y, 8, ttlstr, yd, binc)
    
    ttlstr = 'Estimated states error';
    plot_states(t, x-y, 9, ttlstr, yd, 0);
    
    plot_attitude(t, y, us);
end

function [] = plot_states(t, data, fign, ttlstr, yd, binc)
    figure(fign)
    clf;
    
    subplot(2,2,1);
    bin = 1:3;
    plot(t, data(bin,:));
    grid minor;
    title(sprintf('%s x [m]', ttlstr));
    legend('x','y','z', 'Location', 'best');
    hold on;
    for i = bin
        if(ismember(i, binc))
            h = plot(t, yd(i,:), 'k--');
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
    end
    
    subplot(2,2,2);
    bin = 4:6;
    plot(t, data(bin,:));
    grid minor;
    title(sprintf('%s v [m/s]', ttlstr));
    legend('u','v','w', 'Location', 'best');
    hold on;
    for i = bin
        if(ismember(i, binc))
            h = plot(t, yd(i,:), 'k--');
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
    end
    
    subplot(2,2,3);
    bin = 7:9;
    plot(t, rad2deg(data(bin,:)));
    grid minor;
    title(sprintf('%s %s [deg/s]', ttlstr, '\omega'));
    legend('p','q','r', 'Location', 'best');
    hold on;
    for i = bin
        if(ismember(i, binc))
            h = plot(t, yd(i,:), 'k--');
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
    end
    
    subplot(2,2,4);
    bin = 10:12;
    plot(t, rad2deg(data(bin,:)));
    grid minor;
    title(sprintf('%s %s [deg]', ttlstr, '\epsilon'));
    legend('\phi','\theta','\psi', 'Location', 'best');
    hold on;
    for i = bin
        if(ismember(i, binc))
            h = plot(t, yd(i,:), 'k--');
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
    end    
end

function [] = plot_attitude(ts, ys, us)
    
    global y_pert_e
    Epss = ys(10:12,:)';
    Cs = us';
    
    figure(5000)
    plot(ts, rad2deg(Epss(:,1)), 'b-', 'DisplayName', 'y');
    hold on;
    plot(ts, 0*ts, 'r-', 'DisplayName', 'yN');
    plot(ts, 0*ts + rad2deg(y_pert_e(10))*(1-0.632), 'r--', 'DisplayName', '(100-63.2)%');
    grid minor;
    title('Roll stabilization');
    ylabel('\phi [deg]');
    xlabel('t [s]');
    legend('Location', 'best');
    xlim([4.5, 7.5]);
    axis square;
    ps = get(gcf, 'position');
    ps(4) = 271;
    ps(3) = 287;
    set(gcf, 'position', ps);
    

    figure(5001)
    plot(ts, rad2deg(Epss(:,2)), 'b-', 'DisplayName', 'y');
    hold on;
    plot(ts, 0*ts, 'r-', 'DisplayName', 'yN');
    plot(ts, 0*ts + rad2deg(y_pert_e(11))*(1-0.632), 'r--', 'DisplayName', '(100-63.2)%');
    grid minor;
    title('Pitch stabilization');
    ylabel('\theta [deg]');
    xlabel('t [s]');
    legend('Location', 'best');
    xlim([4.5, 7.5]);
    axis square;
    ps = get(gcf, 'position');
    ps(4) = 271;
    ps(3) = 287;
    set(gcf, 'position', ps);
    

    figure(5002)
    plot(ts, rad2deg(Epss(:,3)), 'b-', 'DisplayName', 'y');
    hold on;
    plot(ts, 0*ts, 'r-', 'DisplayName', 'yN');
    plot(ts, 0*ts + rad2deg(y_pert_e(12))*(1-0.632), 'r--', 'DisplayName', '(100-63.2)%');
    grid minor;
    title('Yaw stabilization');
    ylabel('\psi [deg]');
    xlabel('t [s]');
    legend('Location', 'best');
    xlim([4.5, 7.5]);
    axis square;
    ps = get(gcf, 'position');
    ps(4) = 271;
    ps(3) = 287;
    set(gcf, 'position', ps);
    



    figure(5010)
    plot(ts, Cs(:,1:2), '-');
    hold on;
    plot(ts, Cs(:,3:4), '--');
    grid minor;
    title('Control');
    ylabel('u [ ]');
    xlabel('t [s]');
    legend('u1', 'u2', 'u3', 'u4', 'Location', 'best');
    xlim([4.5, 7.5]);
    axis square;
    ps = get(gcf, 'position');
    ps(4) = 271;
    ps(3) = 287;
    set(gcf, 'position', ps);
    
    
%     path = fullfile('..', '..', 'Report','Imgs2'); 
%     if (y_pert_e(10) == 15) && (y_pert_e(11) == 0) && (y_pert_e(12) == 0)
%         filename = fullfile(path, 'stabilisation-sim_phi');
%         saveas(5000, filename, 'epsc');
%         filename = fullfile(path, 'stabilisation-sim_phi-cntrl');
%         saveas(5010, filename, 'epsc');
%     end
%     if (y_pert_e(10) == 0) && (y_pert_e(11) == 15) && (y_pert_e(12) == 0)
%         filename = fullfile(path, 'stabilisation-sim_theta');
%         saveas(5001, filename, 'epsc');
%         filename = fullfile(path, 'stabilisation-sim_theta-cntrl');
%         saveas(5010, filename, 'epsc');
%     end
%     if (y_pert_e(10) == 0) && (y_pert_e(11) == 0) && (y_pert_e(12) == 15)
%         filename = fullfile(path, 'stabilisation-sim_psi');
%         saveas(5002, filename, 'epsc');
%         filename = fullfile(path, 'stabilisation-sim_psi-cntrl');
%         saveas(5010, filename, 'epsc');
%     end
end


function [] = PX4_printouts(Kht, u0, binx, binu, bine)
    for i = 1:size(Kht,1)
        n = 0;
        for j = binx(1:end)
            fprintf('K_feedback_y(%d,%d) = %8.4ff; ', i-1, n, Kht(i,j));
            n = n + 1;
        end
        fprintf('\n');
    end
    fprintf('\n');

    for i = 1:size(Kht,1)
        n = 0;
        for j = binu
            fprintf('K_feedback_uf(%d,%d) = %8.4ff; ', i-1, n, Kht(i,j));
            n = n + 1;
        end
        fprintf('\n');
    end
    fprintf('\n');

    for i = 1:size(Kht,1)
        n = 0;
        for j = bine
            fprintf('K_feedback_int(%d,%d) = %8.4ff; ', i-1, n, Kht(i,j));
            n = n + 1;
        end
        fprintf('\n');
    end
    fprintf('\n');
    
    for i = 1:length(u0)
        fprintf('u_nominal_control(%d,0) = %5.4ff;\n', i-1, u0(i));
    end
    
end


function [] = plot_eigs(eigs, eigsA, eigsMC, fign, Ts)
    
    figure(fign)
    hold on;
    
    if Ts == 0
        strA = 'A';
        strABK = '(A-BK)';
        strABKMC = '(A-BK) M-C';
    else
        strA = '\Phi';
        strABK = '(\Phi-\GammaK)';
        strABKMC = '(\Phi-\GammaK) M-C';
    end
    
    if ~isempty(eigsMC)
        plot(real(eigsMC), imag(eigsMC), 'b*', 'DisplayName', strABKMC, 'MarkerSize', 3);
    end
    
    if ~isempty(eigsA)
        plot(real(eigsA), imag(eigsA), 'kx', 'DisplayName', strA, 'MarkerSize', 10);
    end
    
    plot(real(eigs), imag(eigs), 'ko', 'DisplayName', strABK, 'MarkerSize', 10);
    
    
    if Ts == 0
        mn = min(real([eigsA; eigs]));
        bin = linspace(mn, -mn, 100);
        h1 = plot(bin, 0*bin, 'r--');
        h2 = plot(0*bin, bin, 'r--');
        set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        title('Eigenvalues');
    else
        ang=0:0.01:2*pi; 
        xp=1*cos(ang);
        yp=1*sin(ang);
        h1 = plot(xp,yp, 'r--');
        set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        axis equal;
        title(sprintf('Eigenvalues Fs=%3.1fHz', 1/Ts));
    end
    grid minor;
    xlabel('real');
    ylabel('imag');
    legend('show');
    set(gcf, 'color', 'white');
    hold off;
end

function [z] = ode1d(Aht, Bht, Kht, z0, t, ud)
    dt = t(2) - t(1);
    sysc = ss(Aht, Bht, 0*Aht, 0*Bht);
    sysd = c2d(sysc, dt);
    z = zeros(length(z0), length(t));
    z(:,1) = z0;
    delay = 0.00; % control delay [s]
    idelay = ceil(delay/dt);
    for i = 1:length(t)-1
        if i-idelay > 0
            u = bound_control(-Kht*z(:,i-idelay), ud(t(i-idelay)));
        else
            u = ud(t(i));
        end
        z(:,i+1) = sysd.A*z(:,i) + sysd.B*u;
    end
end









