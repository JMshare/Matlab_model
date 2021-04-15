%%% Prop thrust model - quadratic %%%

clear all;
close all;
clc;

rho = 1.225; % [kg/m^3]
D = 0.0254*16; % in to [m]

filename = 'my_prop1.txt';
A = importdata(filename);
polars = A.data;
J = polars(:,1);
CT = polars(:,3);
X = fliplr(vander(J));
X = X(:, 1:3);
y = CT;
a = (X'*X)\(X'*y);

cT = @(V,n) max((a(3)*(V./(n*D)).^2 + a(2)*(V./(n*D)) + a(1)), min(CT));
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
plot(J, CT, 'b.');

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
Vs = ((-a(2) - sqrt(a(2)^2 - 4*a(3)*a(1)))/(2*a(3)) * ns*D);
plot(ns*60, Vs, 'r');
grid minor;
xlabel('n [rpm]');
ylabel('Vp [m/s]');
title('Velocity of air passing through propeller');

% omgs = (2*pi)*ns; % to rad/s
% Vs = omgs*D/2;
% hold on;
% plot(ns*60, Vs, 'r--'); % velocity of the propeller tip

figure(5)
A = (pi*D^2)/4;
V = linspace(0, 35, 100);
n = 6000/60; % rpm to [rps]
Vss = @(V,n) sqrt(2*T(V,n)/(rho*A) + V.^2);
plot(V, Vss(V,n), 'r');
grid minor;
hold on;
plot(V, V*0, 'k--');
xlabel('V [m/s]');
ylabel('V slipstream [m/s]');
title('Slipstream velocity(V): n from 0 to 10000 rpm');
ns = [0:1000:10000]/60;
for i=1:length(ns)
    Ts = Vss(V,ns(i));
    plot(V, Ts, '--');
    text(V(1), Ts(1), sprintf(' n=%2.0f',ns(i)*60));
end

figure(6)
n = linspace(0, 10000, 100)/60;
V = 30;
Vs = 0:5:40;
plot(n*60, Vss(V,n), 'r');
grid minor;
hold on;
xlabel('n [rpm]');
ylabel('Vss [m/s]');
title('Slipstream velocity(n): V from 0 to 40 m/s');
for i=1:length(Vs)
    Ts = Vss(Vs(i),n);
    plot(n*60, Ts, '--');
    text(n(end)*60, Ts(end), sprintf('V=%2.0f',Vs(i)));
end


%%% CP power coef
CP = polars(:,4);
X = fliplr(vander(J));
X = X(:, 1:3);
y = CP;
a = (X'*X)\(X'*y);

cP = @(V,n) max((a(3)*(V./(n*D)).^2 + a(2)*(V./(n*D)) + a(1)), min(CP));
P = @(V,n) cP(V,n)*rho.*(n.^3)*D^5;

n = 6000/60; % rpm to [rps]
V = linspace(0, 35, 100);

figure(10)
plot(V/(n*D), cP(V,n), 'b');
grid minor;
hold on;
plot(V/(n*D), V*0, 'k--');
xlabel('J=V/(nD)');
ylabel('cP');
title('Power coef');
plot(J, CP, 'b.');

figure(11)
plot(V, P(V,n), 'r');
grid minor;
hold on;
plot(V, V*0, 'k--');
xlabel('V [m/s]');
ylabel('P [N.m]');
title('Power(V): n from 0 to 10000 rpm');
ns = [0:1000:10000]/60;
for i=1:length(ns)
    Ps = P(V,ns(i));
    plot(V, Ps, '--');
    text(V(1), Ps(1), sprintf(' n=%2.0f',ns(i)*60));
end

figure(12)
n = linspace(0, 10000, 100)/60;
V = 30;
Vs = 0:5:40;
plot(n*60, P(V,n), 'r');
grid minor;
hold on;
xlabel('n [rpm]');
ylabel('P [N.m]');
title('Power(n): V from 0 to 40 m/s');
for i=1:length(Vs)
    Ps = P(Vs(i),n);
    plot(n*60, Ps, '--');
    text(n(end)*60, Ps(end), sprintf('V=%2.0f',Vs(i)));
end


%%% Pressure rise DelP
A = (pi*D^2)/4;
DelP = @(V,n) T(V,n)/A;

n = 6000/60; % rpm to [rps]
V = linspace(0, 35, 100);

figure(20)
plot(V, DelP(V,n), 'r');
grid minor;
hold on;
plot(V, V*0, 'k--');
xlabel('V [m/s]');
ylabel('DelP [Pa]');
title('Pressure-jump(V): n from 0 to 10000 rpm');
ns = [0:1000:10000]/60;
for i=1:length(ns)
    Ps = DelP(V,ns(i));
    plot(V, Ps, '--');
    text(V(1), Ps(1), sprintf(' n=%2.0f',ns(i)*60));
end

figure(21)
n = linspace(0, 10000, 100)/60;
V = 30;
Vs = 0:5:40;
plot(n*60, DelP(V,n), 'r');
grid minor;
hold on;
xlabel('n [rpm]');
ylabel('DelP [Pa]');
title('Pressure-jump(n): V from 0 to 40 m/s');
for i=1:length(Vs)
    Ps = DelP(Vs(i),n);
    plot(n*60, Ps, '--');
    text(n(end)*60, Ps(end), sprintf('V=%2.0f',Vs(i)));
end


% fit the polynomial for fluent
n = 10000/60; % rpm to [rps]
Vmax = 32;
V = linspace(0, Vmax, 100)';
X = fliplr(vander(V));
X = X(:, 1:3);
y = DelP(V,n);
a = (X'*X)\(X'*y);

DelP_fit = @(V) a(3)*V.^2 + a(2)*V + a(1);

figure(22)
plot(V, DelP(V,n), 'b-');
hold on;
plot(V, DelP_fit(V), 'r--');
grid minor;
title('Prssure-jump fit');
xlabel('V [m/s]');
ylabel('DelP [Pa]');
fprintf('DelP = a2.V^2 + a1.V + a0\n');
fprintf('a2 %f, a1 %f, a0 %f\n', a(3), a(2), a(1));


