function [] = Eigenvalues_plot(A, fig, transp)

    global cols
    cols = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]; [0.6350 0.0780 0.1840]];
    
    i_reorder = [1,2,3,4,6,8,11,5,7,9,10,12];
    A_reordered = A(i_reorder, i_reorder);
    A_dyn = A_reordered(4:11, 4:11);
    A_dyn_long = A_dyn(1:4,1:4);
    A_dyn_later = A_dyn(5:8,5:8);
    
    eigs = eig(A_dyn);
    eigs_long = eig(A_dyn_long);
    eigs_later = eig(A_dyn_later);

    
    %% Plot eigenvalues
    figure(fig+1)
    plot_eigs2(eigs, eigs_long, eigs_later, transp);
    
    
 
    
    figure(fig+2)
    plot_eigs_later(A_dyn_later);
    
    figure(fig+3)
    plot_eigs_long(A_dyn_long);
end




function [] = plot_eigs2(eigs, eigs_long, eigs_later, transp)

    global cols
    
    hold on;
    eigs_phug1  = eigs_long(3);
    eigs_phug2  = eigs_long(4);
    eigs_short = eigs_long(1:2);
    eigs_roll = eigs_later([1]);
    eigs_spiral = eigs_later([4]);
    eigs_dutch = eigs_later(2:3);
    h = scatter(real(eigs_phug1), imag(eigs_phug1), 30, cols(1,:), 'filled', 'HandleVisibility', 'off');
    h.MarkerEdgeAlpha = transp;
    h.MarkerFaceAlpha = transp;
    h = scatter(real(eigs_phug2), imag(eigs_phug2), 30, cols(1,:), 'filled', 'HandleVisibility', 'off');
    h.MarkerEdgeAlpha = transp;
    h.MarkerFaceAlpha = transp;
    h = scatter(real(eigs_short), imag(eigs_short), 240, cols(1,:), '*', 'HandleVisibility', 'off', 'LineWidth', 1.5);
    h.MarkerEdgeAlpha = transp;
    h.MarkerFaceAlpha = transp;
    h = scatter(real(eigs_roll), imag(eigs_roll), 240, cols(5,:), 'x', 'HandleVisibility', 'off', 'LineWidth', 2);
    h.MarkerEdgeAlpha = transp;
    h.MarkerFaceAlpha = transp;
    h = scatter(real(eigs_spiral), imag(eigs_spiral), 240, cols(5,:), 'x', 'HandleVisibility', 'off', 'LineWidth', 2);
    h.MarkerEdgeAlpha = transp;
    h.MarkerFaceAlpha = transp;
    h = scatter(real(eigs_dutch), imag(eigs_dutch), 240, cols(5,:), '+', 'HandleVisibility', 'off', 'LineWidth', 1.5);
    h.MarkerEdgeAlpha = transp;
    h.MarkerFaceAlpha = transp;
    mn = -11;
    bin = linspace(mn, -mn, 100);
    plot(bin, 0*bin, 'r--', 'HandleVisibility', 'off');
    plot(0*bin, bin, 'r--', 'HandleVisibility', 'off');
    grid minor;
    xlabel('$\rm{Re}(\lambda) \, [\rm{s}^{-1}]$', 'Interpreter', 'Latex', 'FontSize', 14);
    ylabel('$\rm{Im}(\lambda) \, [\rm{s}^{-1}]$', 'Interpreter', 'Latex', 'FontSize', 14);
    %title('Eigenvalues of the system matrix A');
    legend('show');
    set(gcf, 'color', 'white');
    hold off;
    axis equal;
    xlim([-11, 3]);
    ylim([-11, 11]);
    yticks([-10:5:10]);
    xticks([-10:5:10]);
end




function [] = plot_eigs_later(A_dyn_later)
    eigs = eig(A_dyn_later); % v,p,r,phi
    A_roll = A_dyn_later([2,4], [2,4]);
    A_yaw = A_dyn_later([1,3], [1,3]);
    eigs_roll = eig(A_roll);
    eigs_yaw = eig(A_yaw);

    hold on;
    plot(real(eigs), imag(eigs), 'ko', 'DisplayName', 'total', 'MarkerSize', 10, 'HandleVisibility', 'off');
    plot(real(eigs_roll), imag(eigs_roll), 'b*', 'DisplayName','roll', 'MarkerSize', 10, 'HandleVisibility', 'off');
    plot(real(eigs_yaw), imag(eigs_yaw), 'g*', 'DisplayName','yaw', 'MarkerSize', 10, 'HandleVisibility', 'off');
    mn = min(real([eigs_roll; eigs_yaw; eigs]));
    mn = -10;
    bin = linspace(mn, -mn, 100);
    h1 = plot(bin, 0*bin, 'r--');
    h2 = plot(0*bin, bin, 'r--');
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    grid minor;
    axis equal;
    xlabel('$\rm{Re}(\lambda) \, [\rm{s}^{-1}]$', 'Interpreter', 'Latex', 'FontSize', 14);
    ylabel('$\rm{Im}(\lambda) \, [\rm{s}^{-1}]$', 'Interpreter', 'Latex', 'FontSize', 14);
    %title('Eigenvalues of the system matrix A');
    set(gcf, 'color', 'white');
    hold off;
    xlim([-11, 3]);
    ylim([-11, 11]);
    yticks([-10:5:10]);
    xticks([-10:5:10]);
end

function [] = plot_eigs_long(A_dyn_long)
    eigs = eig(A_dyn_long); % u,w,q,tht  
    A_pitch = A_dyn_long([2,3], [2,3]);
    A_vel = A_dyn_long([1,4], [1,4]);
    eigs_pitch = eig(A_pitch);
    eigs_vel = eig(A_vel);

    hold on;
    plot(real(eigs), imag(eigs), 'ko', 'DisplayName', 'total', 'MarkerSize', 10, 'HandleVisibility', 'off');
    plot(real(eigs_pitch), imag(eigs_pitch), 'b*', 'DisplayName','pitch', 'MarkerSize', 10, 'HandleVisibility', 'off');
    plot(real(eigs_vel), imag(eigs_vel), 'g*', 'DisplayName','vel', 'MarkerSize', 10, 'HandleVisibility', 'off');
    mn = -10;
    bin = linspace(mn, -mn, 100);
    h1 = plot(bin, 0*bin, 'r--');
    h2 = plot(0*bin, bin, 'r--');
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    grid minor;
    axis equal;
    xlabel('$\rm{Re}(\lambda) \, [\rm{s}^{-1}]$', 'Interpreter', 'Latex', 'FontSize', 14);
    ylabel('$\rm{Im}(\lambda) \, [\rm{s}^{-1}]$', 'Interpreter', 'Latex', 'FontSize', 14);
    %title('Eigenvalues of the system matrix A');
    set(gcf, 'color', 'white');
    hold off;
    xlim([-11, 3]);
    ylim([-11, 11]);
    yticks([-10:5:10]);
    xticks([-10:5:10]);
end


