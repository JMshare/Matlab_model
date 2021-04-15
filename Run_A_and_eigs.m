%% Run to get the A matrix and eigenvalue plots.
% Just the results

% clear; % load the Ds from Run_to_trim.m instead
% close all;
% clc;
% format;

%% Load
% Load_params(1);
% Vstar = 20;
% c = [Vstar/2; Vstar/2; 0; 0.5; 0.5];
% [y, u] = Trim_conditions(Vstar, c);
Vstar = 5;
idx = Vstar + 1;
y = Ys(:,idx); u = Us(:,idx);
[V, Alpha] = Extract_aero_states(y);
if idx ==1
    y(1:9) = 0;
end
global c0
c0 = -y(4:6);
u = u + [0;0;0;0;0;0.01;0.01;0.01;0.01]; % can't perturb to -dx if dP at zero already
[A, B, fAero] = ODE_Jac_num(y, u);
A(4:9, 7:12)
fAero(7:9,:)
B(4:9, 6) = - B(4:9, 6) + B(4:9, 7); % left and right motors diff together
B(4:9, 7) = B(4:9, 6)*0.3; % scaled to diff thrust 0.3
Bprnt = B;
Bprnt(7,1) = B(7,1) + B(7,5); % add tailerons to roll
Bprnt([9,7],3) = B([9,7],3) + B([9,7],7); % add diff thrust to yaw
Bprnt(4:9,:)



diag(A(7:9, 7:9))
diag(B(7:9, 1:3))

i_reorder = [1,2,3, 4,6,8,11, 5,7,9,10, 12];
A_reordered = A(i_reorder, i_reorder);
A_dyn = A_reordered(4:11, 4:11);
A_dyn_longit = A_dyn(1:4,1:4);
A_dyn_later = A_dyn(5:8,5:8);
B_reordered = B(i_reorder, :);
B_dyn = B_reordered(4:11,:);

for i = 1:4
    fprintf('%f\t%f\t%f\t%f\n', A_dyn(i,1), A_dyn(i,2), A_dyn(i,3), A_dyn(i,4));
end


%% Plot A matrices with eigenvalues
figure(1)
subplot(1,2,1)
bin = A_reordered;
bin = scale_log(bin);
ttl = sprintf('A matrix reordered');
plot_matrix12(bin, ttl);

subplot(1,2,2)
eigs = eig(A_dyn);
eigs_long = eig(A_dyn_longit);
eigs_later = eig(A_dyn_later);
plot_eigs(eigs, eigs_long, eigs_later);
set(gcf, 'Position', [282   376   747   297]);


%% Plot B matrix
figure(2); clf;
for i = 1:length(Vs)
    y = Ys(:,i); u = Us(:,i);
    global c0
    c0 = -y(4:6);
    u = u + [0;0;0;0;0;0.01;0.01;0.01;0.01]; % can't perturb to -dx if dP at zero already
    [A, B, fAero] = ODE_Jac_num(y, u);
    B(4:9, 6) = - B(4:9, 6) + B(4:9, 7); % left and right motors diff together
    B(4:9, 7) = B(4:9, 6)*1; % scaled to full -1:1 dP range
    B(4:9, 5) = B(4:9, 5)/0.3; % scalled to full -1:1 dE range
    plot_B_matrix(B, i, length(Vs));
end



%% Plot eigs root locus
idxs = [0, 6, 11, 20, 26]+1;
for i = 1:length(Vs)
    idx = i;
    y = Ys(:,idx); u = Us(:,idx);
    if idx == 1
        y(1:9) = 0;
    end
    global c0
    c0 = -y(4:6);
    u = u + [0;0;0;0;0;0.01;0.01;0.01;0.01]; % can't perturb to -dx if dP at zero already
    [A, B, fAero] = ODE_Jac_num(y, u);
    eigs = eig(A);
    i_reorder = [1,2,3,4,6,8,11,5,7,9,10,12];
    A_reordered = A(i_reorder, i_reorder);
    A_dyn = A_reordered(4:11, 4:11);
    A_dyn_longit = A_dyn(1:4,1:4);
    A_dyn_later = A_dyn(5:8,5:8);
    eigs = eig(A_dyn);
    eigs_long = eig(A_dyn_longit);
    eigs_later = eig(A_dyn_later);
    eigs_phug1  = eigs_long(3);
    eigs_phug2  = eigs_long(4);
    eigs_short1 = eigs_long(1);
    eigs_short2 = eigs_long(2);
    eigs_roll = eigs_later([1]);
    eigs_spiral = eigs_later([4]);
    eigs_dutch1 = eigs_later(2);
    eigs_dutch2 = eigs_later(3);
    
    if nnz(i == idxs) > 0
        transp = (i/length(Vs))/1.5+(1-1/1.5);
        Eigenvalues_plot(A, 30, transp);
        
        figure(31); hold on;
        text(real(eigs_short1)-1, imag(eigs_short1)-0.5, string(idx-1), 'HandleVisibility', 'off');
        text(real(eigs_roll), imag(eigs_roll)+0.8, string(idx-1), 'HandleVisibility', 'off');
        text(real(eigs_dutch1)+0.5, imag(eigs_dutch1)+0.5, string(idx-1), 'HandleVisibility', 'off');
        %text(real(eigs_spiral), imag(eigs_spiral)-0.8, string(idx-1), 'HandleVisibility', 'off');
        
        figure(32); hold on;
        text(real(eigs_roll), imag(eigs_roll)+0.8, string(idx-1), 'HandleVisibility', 'off');
        text(real(eigs_dutch1)+0.5, imag(eigs_dutch1)+0.5, string(idx-1), 'HandleVisibility', 'off');
        
        figure(33); hold on;
        text(real(eigs_short1)-1, imag(eigs_short1)-0.5, string(idx-1), 'HandleVisibility', 'off');
        text(real(eigs_phug2), imag(eigs_phug2)+0.8, string(idx-1), 'HandleVisibility', 'off');
    end
    
    figure(31); hold on;
    plot(real(eigs), imag(eigs), 'k.', 'HandleVisibility', 'off');
    
    if i == length(Vs)
        figure(31); hold on;
        scatter([-100],[-100], 30, cols(1,:), 'filled');
        scatter([-100],[-100], 240, cols(1,:), '*');
        scatter([-100],[-100], 240, cols(5,:), '+');
        scatter([-100],[-100], 240, cols(5,:), 'x');
        scatter([-100],[-100], 240, cols(5,:), 'x');
        [~,h] = legend({'Phugoid', 'Short', 'Dutch', 'Roll & Spiral', '# V [m/s]'}, 'FontSize', 11);
        M = findobj(h,'type','patch');
        set(M(2:4),'MarkerSize',sqrt(240));
        set(M(1),'MarkerSize',sqrt(30));
        set(M(5), 'Marker', 'none');
        
        figure(32); hold on;
        plot(-100, -100, 'ko');
        plot(-100, -100, 'b*');
        plot(-100, -100, 'g*');
        legend({'total', 'roll', 'yaw'}, 'FontSize', 10);
        
        figure(33); hold on;
        plot(-100, -100, 'ko');
        plot(-100, -100, 'b*');
        plot(-100, -100, 'g*');
        legend({'total', 'short', 'phugoid'}, 'FontSize', 10);        
    end
end


function plot_B_matrix(B, i, I)
    V = i-1;
    B = B(7:9,:);
    cols = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]; [0.6350 0.0780 0.1840]];
    msz = 20;
    hvis = 'off';
    if i == I
        hvis = 'on';
    end
    subplot(1,3,1)
    hold on;
    plot(V, B(1,1), '.', 'color', cols(1,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    plot(V, B(1,2), '.', 'color', cols(2,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    plot(V, B(1,3), '.', 'color', cols(3,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    plot(V, B(1,5), '.', 'color', cols(4,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    plot(V, B(1,7), '.', 'color', cols(5,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    
    subplot(1,3,2)
    hold on;
    plot(V, B(2,1), '.', 'color', cols(1,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    plot(V, B(2,2), '.', 'color', cols(2,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    plot(V, B(2,3), '.', 'color', cols(3,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    plot(V, B(2,5), '.', 'color', cols(4,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    plot(V, B(2,7), '.', 'color', cols(5,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    
    subplot(1,3,3)
    hold on;
    plot(V, B(3,1), '.', 'color', cols(1,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    plot(V, B(3,2), '.', 'color', cols(2,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    plot(V, B(3,3), '.', 'color', cols(3,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    plot(V, B(3,5), '.', 'color', cols(4,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    plot(V, B(3,7), '.', 'color', cols(5,:), 'HandleVisibility', hvis, 'MarkerSize', msz);
    
    if i == I
        ylm = [-140, 20];
        set(gcf, 'Position', [427   403   722   270]);
        subplot(1,3,1)
        title('Roll moment [Nm]');
        grid minor;
        xlabel('V [m/s]');
        %legend('Ailerons', 'Elevators', 'Rudders', 'Tailerons', 'Diff thrust');
        ylim(ylm);
        xlim([0, 30]);
        subplot(1,3,2)
        title('Pitch moment [Nm]');
        grid minor;
        xlabel('V [m/s]');
        %legend('Ailerons', 'Elevators', 'Rudders', 'Tailerons', 'Diff thrust');
        ylim(ylm);
        xlim([0, 30]);
        subplot(1,3,3)
        title('Yaw moment [Nm]');
        grid minor;
        xlabel('V [m/s]');
        legend({'Ailerons', 'Elevators', 'Rudders', 'Tailerons', 'Diff thrust'}, 'FontSize', 10);
        ylim(ylm);
        xlim([0, 30]);
    end
end



function [] = plot_matrix12(A, ttl)
    mn = min(min(A));
    mx = max(max(A));
    mmx = max([abs(mn), abs(mx)]);
    imagesc(A, 'x', 0.5, 'y', 0.5, [-mmx, mmx]);
    colormap(bluewhitered(256));
    h = colorbar;
    ylabel(h, 'log scaled');
    hold on;
    plot(0:13, 0*(0:13)+3, 'k');
    plot(0:13, 0*(0:13)+7, 'k');
    plot(0:13, 0*(0:13)+11, 'k');
    plot(0*(0:13)+3, 0:13, 'k');
    plot(0*(0:13)+7, 0:13, 'k');
    plot(0*(0:13)+11, 0:13, 'k');
    ax = gca;
    ax.XTick = 0.5:12.5;
    ax.XTickLabel = {'''x','''y','''z','''u','''w','''q','''\theta', '''v', '''p', '''r', '''\phi', '''\psi'};
    ax.YTick = 0.5:12.5;
    ax.YTickLabel = {'dx','dy','dz','du','dw','dq','d\theta', 'dv', 'dp', 'dr', 'd\phi','d\psi'};
    ax.TickLength = [0, 0];
    title(ttl);
    set(gcf, 'color', 'white');
    axis equal;
    ylim([0,13]);
    
end


function [] = plot_eigs(eigs, eigs_long, eigs_later)
    cols = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]; [0.6350 0.0780 0.1840]];
    
    plot(real(eigs), imag(eigs), 'ko', 'DisplayName','total', 'MarkerSize', 10);
    hold on;
    plot(real(eigs_long), imag(eigs_long), '*', 'DisplayName','longitudinal', 'MarkerSize', 10, 'color', cols(1,:));
    plot(real(eigs_later), imag(eigs_later), '*', 'DisplayName','lateral', 'MarkerSize', 10, 'color', cols(5,:));
    mn = -11;
    bin = mn:-mn;
    h1 = plot(bin, 0*bin, 'r--');
    h2 = plot(0*bin, bin, 'r--');
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    grid minor;
    %axis equal;
    %axis([-1,1,-1,1]);
    xlabel('$\rm{Re}(\lambda) \, [\rm{s}^{-1}]$', 'Interpreter', 'Latex', 'FontSize', 14);
    ylabel('$\rm{Im}(\lambda) \, [\rm{s}^{-1}]$', 'Interpreter', 'Latex', 'FontSize', 14);
    title('Eigenvalues of A');
    legend('show');
    set(gcf, 'color', 'white');
    
    axis equal;
    xlim([-11, 3]);
    ylim([-11, 11]);
    yticks([-10:5:10]);
    xticks([-10:5:10]);
    hold off;
end


function [A_sc] = scale_log(A)
    A_sc = log(abs(A)+1).*sign(A);
end




