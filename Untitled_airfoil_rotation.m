clear all;
close all;
format;


XY = [
    1.0000     0.0013
    0.9500     0.0147
    0.9000     0.0271
    0.8000     0.0489
    0.7000     0.0669
    0.6000     0.0814
    0.5000     0.0919
    0.4000     0.0980
    0.3000     0.0976
    0.2500     0.0941
    0.2000     0.0880
    0.1500     0.0789
    0.1000     0.0659
    0.0750     0.0576
    0.0500     0.0473
    0.0250     0.0339
    0.0125     0.0244
    0.0000     0.0000
    0.0125     -0.0143
    0.0250     -0.0195
    0.0500     -0.0249
    0.0750     -0.0274
    0.1000     -0.0286
    0.1500     -0.0288
    0.2000     -0.0274
    0.2500     -0.0250
    0.3000     -0.0226
    0.4000     -0.0180
    0.5000     -0.0140
    0.6000     -0.0100
    0.7000     -0.0065
    0.8000     -0.0039
    0.9000     -0.0022
    0.9500     -0.0016
    1.0000     -0.0013];


figure(1)
grid minor;
hold on;
axis equal;
    
L_proj = plot_airfoil_projection(XY);


alpha_min = fmincon(@(alpha) rotate_airfoil(alpha, XY), 0);

alpha = alpha_min;
R = [cosd(alpha), -sind(alpha); sind(alpha), cosd(alpha)];
XY_rot = (R*XY')';
L_proj_rot = plot_airfoil_projection(XY_rot);
[L_proj, L_proj_rot]


function [L_proj, XY] = rotate_airfoil(alpha, XY)
    R = [cosd(alpha), -sind(alpha); sind(alpha), cosd(alpha)];
    XY = (R*XY')';
    [mn, imn] = min(XY(:,2));
    [mx, imx] = max(XY(:,2));
    L_proj = mx - mn;
end

function [L_proj] = plot_airfoil_projection(XY)
    [mn, imn] = min(XY(:,2));
    [mx, imx] = max(XY(:,2));
    L_proj = mx - mn;

    plot(XY(:,1), XY(:,2));
    
    plot([XY(imn,1), XY(imn,1)], [mn, mn+L_proj], 'r--');
    plot([XY(imn,1), XY(imx,1)], [mx, mx], 'r--')
end



