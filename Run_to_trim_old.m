%%% Run trimming %%%


clear;
close all;
clc;
format;


Load_params(1);
fign = 1000;
Vs = 0:1:30;
Vs(1) = 0.001;
cslast = [Vs(1)/2; Vs(1)/2; 0; 0.5; 0.5];
Us = [];
Ys = [];
for i = 1:length(Vs)
    v = Vs(i);
    [y, u, cslast] = Trim_conditions(v, cslast);
    Us = [Us;u'];
    Ys = [Ys;y'];
    [F, M, hp, Fa, Fg, Ft, Fcc, Fcc_stall, Ma, Mt, Mcc, Mcc_stall] = Forces_and_Moments(y,u);
    [V, Alpha] = Extract_aero_states(y);
    
    [La, Da, Ma] = To_wind_frame(Fa, Ma, Alpha);
    [Lg, Dg, ~] = To_wind_frame(Fg, [0;0;0], Alpha);
    [Lt, Dt, Mt] = To_wind_frame(Ft, Mt, Alpha);
    [Lcc, Dcc, Mcc] = To_wind_frame(Fcc, Mcc, Alpha);
    [Lcc_s, Dcc_s, Mcc_s] = To_wind_frame(Fcc_stall, Mcc_stall, Alpha);
    
    figure(fign)
    msz = 20;
    subplot(1,5,1);
    hold on;
    plot(v, La, 'b.', 'MarkerSize', msz);
    plot(v, Lg, 'g.', 'MarkerSize', msz);
    plot(v, Lt, 'r.', 'MarkerSize', msz);
    plot(v, Lcc, 'k.', 'MarkerSize', msz);
    plot(v, Lcc_s, 'k*');
    if abs(La+Lg+Lt+Lcc+Lcc_s) > 1
        plot(v, La+Lg+Lt+Lcc+Lcc_s, 'y*', 'MarkerSize', msz);
    end
    subplot(1,5,2);
    hold on;
    plot(v, Da, 'b.', 'MarkerSize', msz);
    plot(v, Dt, 'r.','MarkerSize', msz);
    plot(v, Dcc, 'k.', 'MarkerSize', msz);
    plot(v, Dcc_s, 'k*');
    if abs(Da+Dt+Dcc+Dcc_s) > 1
        plot(v, Da+Dt+Dcc+Dcc_s, 'y*', 'MarkerSize', msz);
    end
    
    subplot(1,5,3);
    hold on;
    plot(v, Ma, 'b.','MarkerSize', msz);
    plot(v, Mt, 'r.', 'MarkerSize', msz);
    plot(v, Mcc, 'k.', 'MarkerSize', msz);
    plot(v, Mcc_s, 'k*');
    if abs(Ma+Mt+Mcc+Mcc_s) > 1
        plot(v, Ma+Mt+Mcc+Mcc_s, 'y*', 'MarkerSize', msz);
    end
    
    subplot(1,5,4);
    hold on;
    plot(v, Alpha, 'b.','MarkerSize', msz);
    
    subplot(1,5,5);
    hold on;
%     global P1 P2 P3 U_max dP1min dP2min dP3min
%     power_front = norm(P3(v,u(8)*U_max(8,8))*(u(8)>dP3min(v)))/norm(P3(v,U_max(8,8)));
%     power_ducts = (norm(P1(v,u(6)*U_max(6,6))*(u(6)>dP1min(v))) ...
%                  + norm(P2(v,u(7)*U_max(7,7))*(u(7)>dP2min(v)))) ...
%                  /(norm(P1(v,U_max(6,6))) + norm(P2(v,U_max(7,7))));
%     power_total =  (power_front + power_ducts)/3;
%     plot(v, power_ducts, 'k.', 'MarkerSize', msz);
%     plot(v, power_front, 'r.', 'MarkerSize', msz);
    plot(v, (u(6)+u(7))/2, 'k.', 'MarkerSize', msz);
    plot(v, u(8), 'r.', 'MarkerSize', msz);
    plot(v, u(2), 'b.', 'MarkerSize', msz);
%     plot(v, u(4), 'g.', 'MarkerSize', msz);
    
    if norm([F(1); F(3); M(2)]) < 1
        figure(fign+1)
        subplot(1,2,1);
        hold on;
        plot(v, Alpha, 'b.','MarkerSize', msz);
        subplot(1,2,2);
        hold on;
        plot(v, (u(6)+u(7))/2, 'k.', 'MarkerSize', msz);
        plot(v, u(8), 'r.', 'MarkerSize', msz);
        plot(v, u(2), 'b.', 'MarkerSize', msz);
%         plot(v, u(4), 'g.', 'MarkerSize', msz);

        figure(fign+2)
        subplot(1,3,1);
        hold on;
        plot(v, La, 'b.', 'MarkerSize', msz);
        plot(v, Lg, 'g.', 'MarkerSize', msz);
        plot(v, Lt, 'r.', 'MarkerSize', msz);
        plot(v, Lcc, 'k.', 'MarkerSize', msz);
        plot(v, Lcc_s, 'k*');
        subplot(1,3,2);
        hold on;
        plot(v, Da, 'b.', 'MarkerSize', msz);
        plot(v, Dt, 'r.','MarkerSize', msz);
        plot(v, Dcc, 'k.', 'MarkerSize', msz);
        plot(v, Dcc_s, 'k*');
        subplot(1,3,3);
        hold on;
        plot(v, Ma, 'b.','MarkerSize', msz);
        plot(v, Mt, 'r.', 'MarkerSize', msz);
        plot(v, Mcc, 'k.', 'MarkerSize', msz);
        plot(v, Mcc_s, 'k*');
    end
    
end
figure(fign)
subplot(1,5,1);
xlabel('V [m/s]');
set(gca, 'xdir', 'reverse');
grid minor;
title('L');
legend('La', 'Lg', 'Lt', 'Lcc', 'Lcc_s', 'location', 'best');
subplot(1,5,2);
xlabel('V [m/s]');
set(gca, 'xdir', 'reverse');
grid minor;
title('D');
legend('Da', 'Dt', 'Dcc', 'Dcc_s', 'location', 'best');
subplot(1,5,3);
xlabel('V [m/s]');
set(gca, 'xdir', 'reverse');
grid minor;
title('M');
legend('Ma', 'Mt', 'Mcc', 'Mcc_s','location', 'best');
if norm(get(gca, 'Ylim')) < 10
    ylim([-10, 10]);
end
subplot(1,5,4);
xlabel('V [m/s]');
set(gca, 'xdir', 'reverse');
grid minor;
title('Alpha');
subplot(1,5,5);
xlabel('V [m/s]');
set(gca, 'xdir', 'reverse');
grid minor;
title('Controls');
legend('u_{PW}', 'u_{PF}', 'u_E', 'u_{AS}', 'location', 'best');


figure(fign+1)
set(gcf, 'position', [750, 376, 664, 370/1.5]);
subplot(1,2,1);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('Alpha');
figure(fign+1)
subplot(1,2,2);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('Controls');
legend('u_{PW}', 'u_{PF}', 'u_E', 'u_{AS}', 'location', 'best');

figure(fign+2)
set(gcf, 'position', [750, 376, 664, 370/1.5]);
subplot(1,3,1);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('Lift L');
legend('La', 'Lg', 'Lt', 'Lcc', 'Lcc_s', 'location', 'best');
subplot(1,3,2);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('Drag D');
legend('Da', 'Dt', 'Dcc', 'Dcc_s', 'location', 'best');
subplot(1,3,3);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('Pitch m');
legend('ma', 'mt', 'mcc', 'mcc_s','location', 'best');
if norm(get(gca, 'Ylim')) < 10
    ylim([-10, 10]);
end



function [L, D, M] = To_wind_frame(F, M, Alpha)
    [R_V] = Wind_frame_matrix(Alpha);
    FM = [F;M];
    FMw = R_V*FM;
    L = FMw(1);
    D = FMw(3);
    M = M(2);
end


