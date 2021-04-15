%%% Prop thrust model - linear %%%

clear all;
close all;
clc;

a1 = -0.05/0.3;
a0 = 0.125;
mincT = -0.02;
rho = 1.225; % [kg/m^3]
D = 0.0254*16; % in to [m]

cT = @(V,n) max(a1*(V./(n*D)) + a0, mincT);
T = @(V,n) cT(V,n)*rho.*(n.^2)*D^4;
%T = @(V,n) (a1*V + a0.*n*D).*n*D^3*rho;

n = 6000/60; % rpm to [rps]
V = linspace(0, 35, 100);

figure(1)
plot(V/(n*D), cT(V,n), 'b');
grid minor;
hold on;
plot(V/(n*D), V*0, 'k--');
xlabel('J=V/(nD)');
ylabel('cT');
title('Thrust coef');

figure(2)
plot(V, T(V,n), 'r');
grid minor;
hold on;
plot(V, V*0, 'k--');
xlabel('V [m/s]');
ylabel('T [N]');
title('Thrust(V): n from 0 to 10000 rpm');
ns = [0:1000:10000]/60;
for i=1:length(ns)
    Ts = T(V,ns(i));
    plot(V, Ts, '--');
    text(V(1), Ts(1), sprintf(' n=%2.0f',ns(i)*60));
end

figure(3)
n = linspace(0, 10000, 100)/60;
V = 30;
Vs = 0:5:40;
plot(n*60, T(V,n), 'r');
grid minor;
hold on;
xlabel('n [rpm]');
ylabel('T [N]');
title('Thrust(n): V from 0 to 40 m/s');
for i=1:length(Vs)
    Ts = T(Vs(i),n);
    plot(n*60, Ts, '--');
    text(n(end)*60, Ts(end), sprintf('V=%2.0f',Vs(i)));
end

figure(4)
ns = [0:1000:10000]/60;
Vs = -a0*ns*D/a1;
plot(ns*60, Vs, 'r');
grid minor;
xlabel('n [rpm]');
ylabel('Vp [m/s]');
title('Velocity of air passing through propeller');

% omgs = (2*pi)*ns; % to rad/s
% Vs = omgs*D/2;
% hold on;
% plot(ns*60, Vs, 'r--'); % velocity of the propeller tip





