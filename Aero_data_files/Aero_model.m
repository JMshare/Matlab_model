%%% Newtonian flow aero model approximation %%%

clear;
close all;



%% functions
c = 0.3367; % reference chord [m]
filename = 'T1-30_0 m_s-VLM2-15_0kg-x39_1mm-0.txt';
delimiterIn = ' ';
headerlinesIn = 8;
A = importdata(filename,delimiterIn,headerlinesIn);
polars = A.data;
Alpha = polars(:,1);
XCP = interp1(Alpha, polars(:,end), 0);
polars = polars(:, [3,6,7,8,9,10]);

CLaoa = @(a) 2*sin(a).^2.*cos(a);
CDaoa = @(a) 2*sin(a).^3;
CMaoa = @(a) -2*sin(a).^2*abs(XCP/c);
CL = @(A) interp1(Alpha, polars(:,1), A);
CD = @(A) interp1(Alpha, polars(:,2), A);
CM = @(A) interp1(Alpha, polars(:,5), A);
figure(400)
subplot(1,3,1);
A = linspace(-4, 15, 50);
plot(A, CL(A), 'b');
hold on;
Aaoa = linspace(-0, 90, 50);
plot(Aaoa, CLaoa(deg2rad(Aaoa)), 'b--');
grid minor;
title('CL VLM vs Newtonian');
xlabel('Alpha [deg]');
ylabel('CL');
subplot(1,3,2);
plot(A, CD(A), 'b');
hold on;
plot(Aaoa, CDaoa(deg2rad(Aaoa)), 'b--');
grid minor;
title('CD VLM vs Newtonian');
xlabel('Alpha [deg]');
ylabel('CD');
subplot(1,3,3);
plot(A, CM(A), 'b');
hold on;
plot(Aaoa, CMaoa(deg2rad(Aaoa)), 'b--');
grid minor;
title('CM VLM vs Newtonian');
xlabel('Alpha [deg]');
ylabel('CM');
