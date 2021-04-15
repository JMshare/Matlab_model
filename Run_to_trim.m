%%% Run trimming %%%


% clear;
% close all;
% clc;
% format;

global U_max P1 P3
Load_params(1);
fign = 1000;
Vs = 0:1:30;
Vs(1) = 0.001;
cslast = [Vs(1)/2; Vs(1)/2; 0; 0.5; 0.5];
Us = [];
Ys = [];
Ds = [];
for i = 1:length(Vs)
    v = Vs(i);
    [y, u, cslast] = Trim_conditions(v, cslast);
    Us = [Us,u];
    Ys = [Ys,y];
    [F, M, hp, Fa, Fg, Ft_components, Fcc, Fcc_stall, Ma, Mt_components, Mcc, Mcc_stall, Fa_tail, Ma_tail, Vss, Vssf] = Forces_and_Moments(y,u);
    [V, Alpha] = Extract_aero_states(y);
    
    [La, Da, Ma] = To_wind_frame(Fa, Ma, Alpha);
    [Lg, Dg, ~] = To_wind_frame(Fg, [0;0;0], Alpha);
    Ft_front = Ft_components(:,3);
    Ft_wings = Ft_components(:,1) + Ft_components(:,2);
    Mt_front = Mt_components(:,3);
    Mt_wings = Mt_components(:,1) + Mt_components(:,2);
    [Lt_front, Dt_front, Mt_front] = To_wind_frame(Ft_front, Mt_front, Alpha);
    [Lt_wings, Dt_wings, Mt_wings] = To_wind_frame(Ft_wings, Mt_wings, Alpha);
    [Lcc, Dcc, Mcc] = To_wind_frame(Fcc, Mcc, Alpha);
    [Lcc_s, Dcc_s, Mcc_s] = To_wind_frame(Fcc_stall, Mcc_stall, Alpha);
    [La_tail, Da_tail, Ma_tail] = To_wind_frame(Fa_tail, Ma_tail, Alpha);
    
    dP1 = u(6)*U_max(6,6);
    dP3 = u(8)*U_max(8,8);
    p1 = P1(V, dP1);
    p3 = P3(V, dP3);
    p1max = 1250;
    p3max = 2060;
    pf1 = p1/p1max;
    pf3 = p3/p3max;
    Ds = [Ds; [V, Alpha, La, Da, Lt_front, Dt_front, Lt_wings+Lcc+Lcc_s, Dt_wings+Dcc+Dcc_s, u(2)*U_max(2,2), La_tail, u(6), u(8), Vss, Vssf, pf1, pf3]];
    
    figure(fign)
    msz = 20;
    subplot(1,5,1);
    hold on;
    plot(v, La, 'b.', 'MarkerSize', msz);
    plot(v, Lg, 'g.', 'MarkerSize', msz);
    plot(v, Lt_front, 'r.', 'MarkerSize', msz);
    plot(v, Lt_wings, 'k.', 'MarkerSize', msz);
    plot(v, Lcc_s+Lcc, 'k*');
    if abs(La+Lg+Lt_front+Lt_wings+Lcc+Lcc_s) > 1
        plot(v, La+Lg+Lt_front+Lt_wings+Lcc+Lcc_s, 'y*', 'MarkerSize', msz);
    end
    subplot(1,5,2);
    hold on;
    plot(v, Da, 'b.', 'MarkerSize', msz);
    plot(v, Dt_front, 'r.','MarkerSize', msz);
    plot(v, Dt_wings, 'k.', 'MarkerSize', msz);
    plot(v, Dcc_s + Dcc, 'k*');
    if abs(Da+Dt_front+Dt_wings+Dcc+Dcc_s) > 1
        plot(v, Da+Dt_front+Dt_wings+Dcc+Dcc_s, 'y*', 'MarkerSize', msz);
    end
    
    subplot(1,5,3);
    hold on;
    plot(v, Ma, 'b.','MarkerSize', msz);
    plot(v, Mt_front, 'r.', 'MarkerSize', msz);
    plot(v, Mt_wings, 'k.', 'MarkerSize', msz);
    plot(v, Mcc_s + Mcc, 'k*');
    if abs(Ma+Mt_front+Mt_wings+Mcc+Mcc_s) > 1
        plot(v, Ma+Mt_front+Mt_wings+Mcc+Mcc_s, 'y*', 'MarkerSize', msz);
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
        plot(v, Lt_front, 'r.', 'MarkerSize', msz);
        plot(v, Lt_wings, 'k.', 'MarkerSize', msz);
        plot(v, Lcc_s + Lcc, 'k*');
        subplot(1,3,2);
        hold on;
        plot(v, Da, 'b.', 'MarkerSize', msz);
        plot(v, Dt_front, 'r.','MarkerSize', msz);
        plot(v, Dt_wings, 'k.', 'MarkerSize', msz);
        plot(v, Dcc_s + Dcc, 'k*');
        subplot(1,3,3);
        hold on;
        plot(v, Ma, 'b.','MarkerSize', msz);
        plot(v, Mt_front, 'r.', 'MarkerSize', msz);
        plot(v, Mt_wings, 'k.', 'MarkerSize', msz);
        plot(v, Mcc_s + Mcc, 'k*');
        
        figure(fign+3)
        subplot(1,2,1);
        hold on;
        plot(v, La, 'b.', 'MarkerSize', msz);
        plot(v, Lg, 'g.', 'MarkerSize', msz);
        plot(v, Lt_front, 'r.', 'MarkerSize', msz);
        plot(v, Lt_wings, 'k.', 'MarkerSize', msz);
        plot(v, Lcc_s + Lcc, 'k*');
        subplot(1,2,2);
        hold on;
        plot(v, Da, 'b.', 'MarkerSize', msz);
        plot(v, Dt_front, 'r.','MarkerSize', msz);
        plot(v, Dt_wings, 'k.', 'MarkerSize', msz);
        plot(v, Dcc_s + Dcc, 'k*');
    end
    
end
figure(fign)
subplot(1,5,1);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('L');
legend('La', 'Lg', 'Lt-front', 'Lt-wings', 'Lcc', 'location', 'best');
subplot(1,5,2);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('D');
legend('Da', 'Dt-front', 'Dt-wings', 'Dcc', 'location', 'best');
subplot(1,5,3);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('M');
legend('Ma', 'Mt-front', 'Mt-wings', 'Mcc','location', 'best');
if norm(get(gca, 'Ylim')) < 10
    ylim([-10, 10]);
end
subplot(1,5,4);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('Alpha');
subplot(1,5,5);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('Controls');
legend('u_{PW}', 'u_{PF}', 'u_E', 'u_{AS}', 'location', 'best');


figure(fign+1)
set(gcf, 'position', [750, 376, 664, 370/1.5]);
subplot(1,2,1);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('Alpha [deg]');
figure(fign+1)
subplot(1,2,2);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('Controls [-]');
legend('Wing motors', 'Front engine', 'Elevators', 'location', 'best');

figure(fign+2)
set(gcf, 'position', [750, 376, 664, 370/1.5]);
subplot(1,3,1);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('Lift L [N]');
legend('Aerodynamic', 'Gravitational', 'Front engine', 'Wing motors', 'Channel wings', 'location', 'best');
subplot(1,3,2);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('Drag D [N]');
%legend('Da', 'Dp-front', 'Dp-wings', 'Dcc', 'location', 'best');
subplot(1,3,3);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('Pitch m [Nm]');
% legend('ma', 'mp-front', 'mp-wings', 'mcc','location', 'best');
if norm(get(gca, 'Ylim')) < 10
    ylim([-10, 10]);
end


figure(fign+3)
set(gcf, 'position', [750, 376, 664, 370/1.5]);
subplot(1,2,1);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('Lift L [N]');
legend('Aerodynamic', 'Gravitational', 'Front engine', 'Wing motors', 'Channel wings', 'location', 'best');
subplot(1,2,2);
xlabel('V [m/s]');
%set(gca, 'xdir', 'reverse');
grid minor;
title('Drag D [N]');
% legend('Da', 'Dp-front', 'Dp-wings', 'Dcc', 'location', 'best');



Us2 = [];
Ds2 = [];
Hs2 = [];
Fs2 = [];
Ms2 = [];
Mas2 = [];
Mtfs2 = [];
Mtes2 = [];
for i = 1:length(Vs)
    v = Vs(i);
    u = Us(:,i);
    y = Ys(:,i);
    [u2] = Trim_conditions_later(y, u);
    Us2 = [Us2,u2];
    [F, M, hp, Fa, Fg, Ft_components, Fcc, Fcc_stall, Ma, Mt_components, Mcc, Mcc_stall, Fa_tail, Ma_tail, Vss, Vssf] = Forces_and_Moments(y,u2);
    Mt_front = Mt_components(:,3);
    Mt_wings = Mt_components(:,1) + Mt_components(:,2);
    
    Hs2 = [Hs2,hp];
    Fs2 = [Fs2,F];
    Ms2 = [Ms2,M];
    Mas2 = [Mas2,Ma];
    Mtfs2 = [Mtfs2,Mt_front];
    Mtes2 = [Mtes2,Mt_wings];
    
end
if(abs(Us2(3,2)) > 0.99)
    Us2(3,1) = 1*sign(Us2(3,2)); % rudders completely ineffective for the trim routine to pick it up but we can say it is trying to compensate
end
if(abs(Us2(1,2)) > 0.99)
    Us2(1,1) = 1*sign(Us2(1,2)); % rudders completely ineffective for the trim routine to pick it up but we can say it is trying to compensate
end

is = 1:length(Vs);
is = is-1;
cols = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]; [0.6350 0.0780 0.1840]];
figure()
subplot(1,3,1)
hold on;
plot(is, Mas2(1,:), 'b.','MarkerSize', msz);
plot(is, Mtfs2(1,:), 'r.','MarkerSize', msz);
plot(is, Mtes2(1,:), 'k.','MarkerSize', msz);
title('Roll moment [Nm]');
grid minor;
xlabel('V [m/s]');
ylim([-3,3]);
subplot(1,3,2)
hold on;
plot(is, Mas2(3,:), 'b.','MarkerSize', msz);
plot(is, Mtfs2(3,:), 'r.','MarkerSize', msz);
plot(is, Mtes2(3,:), 'k.','MarkerSize', msz);
title('Yaw moment [Nm]');
grid minor;
xlabel('V [m/s]');
legend('Aerodynmic', 'Front engine', 'Wing motors');
ylim([-3, 3]);
subplot(1,3,3)
hold on;
plot(is, Us2(1,:), '.','MarkerSize', msz, 'color', cols(1,:));
plot(is, Us2(3,:), '.','MarkerSize', msz, 'color', cols(3,:));
title('Controls [-]');
grid minor;
xlabel('V [m/s]');
legend('Ailerons', 'Rudders');
ylim([-1,1]);
set(gcf, 'Position', [427   403   722   270]);


return
%% scale the f_cw effect if you like
is = 1:length(Vs);
bin = is <= 15;
figure()
plot(is(bin)-1, Us0(6,bin), 'ko', 'MarkerSize', msz/2);
hold on;
plot(is(bin)-1, Us(6,bin), 'b.', 'MarkerSize', msz);
grid minor;
title('Wing motors throttle [-]');
legend({'original', 'improved'}, 'FontSize', 11);

yticks([-1, -0.5, 0, 0.5, 1]);
ylim([-0.5, 1]);

Alphas = [];
Alphas0 = [];
for i = 1:length(Vs)
    y = Ys(:,i);
    [V, Alpha] = Extract_aero_states(y);
    Alphas = [Alphas, Alpha];
    
    y0 = Ys0(:,i);
    [V, Alpha0] = Extract_aero_states(y0);
    Alphas0 = [Alphas0, Alpha0];
end
is = 1:length(Vs);
bin = is <= 15;
figure()
plot(is(bin)-1, Alphas0(bin), 'ko', 'MarkerSize', msz/2);
hold on;
plot(is(bin)-1, Alphas(bin), 'b.', 'MarkerSize', msz);
grid minor;
title('Alpha [deg]');
legend({'original', 'improved'}, 'FontSize', 11);


