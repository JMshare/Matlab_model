function Load_params_old(cs)
%% Specify all the parameters of the system here.
% Allows for different case studies selected by cs.
% Parameters are loaded globally after call.
    
    global m I minv Iinv S c b XCP Vcruise
    global Ip1 Ip2 Ip3 Ip4 xp1 xp2 xp3 xp4
    global g rho
    global U_max
    
    %set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on')
    set(groot, 'defaultFigureColor', 'white');
    
    switch cs
        case 1
            %% Mass and geometry
            m = 16.0; % mass [kg]
            Ix = 4.03; % [kg.m^2]
            Iy = 5.57; % [kg.m^2]
            Iz = 9.25; % [kg.m^2]
            Ixz = -0.33; % [kg.m^2]
            I = [Ix   0   -Ixz;
                 0    Iy  0;
                 -Ixz 0   Iz];
            S = 1; % wing area [m^2]
            b = 3; % wing span [m]
            c = 0.35; % MAC [m]
            
            %% Atmo
            rho = 1.225; % kg/m^3]
            g = 9.81; % [m/s^2]
            
            %% Aerodynamics
            Vcruise = 25; % cruise speed [m/s]
            % static and stability
            XCP = -0.35/4;
            filenamestatic.CL = 'Aero_data_files3\CL.txt';
            filenamestatic.CD = 'Aero_data_files3\CD.txt';
            filenamea.SD = 'Aero_data_files3\T7-raw-out0.txt';
            filenamea.nominal = 'Aero_data_files3\T1-25_0 m_s-VLM2-15_0kg-x39_1mm.txt';
            % control
            udE_max = 30; % [deg] max deflection of control surfaces
            udA_max = 30;
            udR_max = 30;
            udES_max = 15;
            udAS_max = 20;
            filenamea.E = 'Aero_data_files3\T1-25_0 m_s-VLM2-15_0kg-x39_1mm-dE03.txt';
            filenamea.R = 'Aero_data_files3\T1-25_0 m_s-VLM2-15_0kg-x39_1mm-dR03.txt';
            filenamea.A = 'Aero_data_files3\T1-25_0 m_s-VLM2-15_0kg-x39_1mm-dA03.txt';
            filenamea.AS = 'Aero_data_files3\T1-25_0 m_s-VLM2-15_0kg-x39_1mm-dAS03.txt';
            filenamea.ES = 'Aero_data_files3\T1-25_0 m_s-VLM2-15_0kg-x39_1mm-dES03.txt';
            
            %% Custer effect
            filenamecc.lift = 'Aero_data_files3\custer_effect.txt';
            filenamecc.stall.CL = 'Aero_data_files3\CLCuster.txt';
            filenamecc.stall.CD = 'Aero_data_files3\CDCuster.txt';
            
            %% Propellers
            filenamep.p1 = 'Aero_data_files3\my_prop_wing.txt';
            filenamep.p2 = 'Aero_data_files3\my_prop_wing.txt';
            filenamep.p3 = 'Aero_data_files3\my_prop_front.txt';
            filenamep.p4 = 'Aero_data_files3\my_prop_front.txt'; % just for plotting
            
            mp1 = 0.072; % [kg] % left
            Rp1 = 0.356/2; % [m] % 0.34 the shorter one but in WT it made no difference
            ep1 = [0,deg2rad(8),0]; % [rad]
            xp1 = [0*-0.35,-0.332,0]; % [m]
            np1_max = 9000/60; % [rps]
            
            mp2 = 0.072; % right
            Rp2 = 0.356/2;
            ep2 = [0,deg2rad(8),0];
            xp2 = [0*-0.35,0.332,0];
            np2_max = 9000/60;
            
            mp3 = 0.09; % front
            Rp3 = 0.457 /2;
            ep3 = [0,0,0];
            xp3 = [0.8,0,0];
            np3_max = 8000/60;
            
            mp4 = 0; % for plotting
            Rp4 = 0.17;
            ep4 = [0,deg2rad(0),0];
            xp4 = [0,0,0];
            np4_max = 9000/60;
        
        otherwise
            fprintf('Error: Case %d not defined.\n', cs);
    end
    
    
    %%% Processing
    % Read the aerodynamic derivatives data
    load_aerodynamic_derivatives(filenamea);
    
    % Read the static aerodynamic coefs with the custer stall effects
    Aero_static_fits(filenamestatic, filenamecc);
    
    
    minv = 1/m;
    Iinv = inv(I);
    
    % Propeller inertial terms
    [Ip1] = Mom_Inertia_estim_prop(mp1, Rp1, ep1);
    [Ip2] = Mom_Inertia_estim_prop(mp2, Rp2, ep2);
    [Ip3] = Mom_Inertia_estim_prop(mp3, Rp3, ep3);
    [Ip4] = Mom_Inertia_estim_prop(mp4, Rp4, ep4);
    
    % Slipstream velocity decay factor
    load_slipstream_decay_factor();
    
    % Propeller thrust data
    load_propeller_thrust_model(Rp1, Rp2, Rp3, Rp4, ep1, ep2, ep3, ep4, np1_max, np2_max, np3_max, np4_max, filenamep);
    
    % Custer effects
    load_Custer_Channel_lift_effect(filenamecc);
    
    % Limits on the controls
    U_max = diag([udA_max; udE_max; udR_max; udAS_max; udES_max; np1_max; np2_max; np3_max; np4_max]);
    
    %%% Check the model
    plot_aero_derivatives();
    plot_thrust_model();
    plot_decay_factor();
    plot_Custer_Channel_lift_effect(filenamecc);
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function definitions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ip] = Mom_Inertia_estim_prop(m, R, ep)
    Ix = 2*(m*R^2)/3;
    Ip = zeros(3,3);
    Ip(1,1) = Ix;
    [~,~,~,~, R, ~,~,~] = Rotation_and_Euler_Matrices(ep);
    Ip = R*Ip;
end


function [] = load_aerodynamic_derivatives(filenamea)
    
    global CLq Cmq CYb CYr Clp Cnb Cnr
    global CYp Clb Clr Cnp
    global CLdE CDdE CmdE CYdR CldR CndR CldA CLdAS CDdAS CmdAS CldES
        
    %% Aerodynamic Derivatives 
    delimiterIn = ',';
    headerlinesIn = 1;
    A = importdata(filenamea.SD,delimiterIn,headerlinesIn);
    polars = A.data;
    Alpha0 = 0;
    polars = interp1(polars(:,1), polars, Alpha0);
    CLq0 = polars(31); 
    Cmq0 = polars(34); 
    CYb0 = polars(35); 
    CYp0 = 0*polars(36); % no effect
    CYr0 = polars(37); 
    Clb0 = 0*polars(38); % no effect
    Clp0 = polars(39); 
    Clr0 = 0*polars(40); % no effect
    Cnb0 = polars(41); 
    Cnp0 = 0*polars(42); % no effect
    Cnr0 = polars(43); 
    
    if not(isempty(strfind(filenamea.SD, 'Aero_data_files3')))
        Cnr0 = 1*Cnr0; % you can modify a coefs value here for case2 comparison
    end

    %% Control Derivatives
    [polars, Alpha, alpha, txt] = read_data_static(filenamea.nominal);
    [polarsE, AlphaE, alphaE, txtE] = read_data_static(filenamea.E);
    [Al, al, polars, polarsE] = interpolate_alpha(Alpha, AlphaE, polars, polarsE, Alpha0);
    polars_dE = (polarsE - polars)/(3);
    CLdE0 = polars_dE(1);
    CmdE0 = polars_dE(5);
    CDdE0 = polars_dE(2);

    [polars, Alpha, alpha, txt] = read_data_static(filenamea.nominal);
    [polarsR, AlphaR, alphaR, txtR] = read_data_static(filenamea.R);
    [Al, al, polars, polarsR] = interpolate_alpha(Alpha, AlphaR, polars, polarsR, Alpha0);
    polars_dR = (polarsR - polars)/(3);
    CYdR0 = polars_dR(3);
    CldR0 = polars_dR(4);
    CndR0 = polars_dR(6);

    [polars, Alpha, alpha, txt] = read_data_static(filenamea.nominal);
    [polarsA, AlphaA, alphaA, txtA] = read_data_static(filenamea.A);
    [Al, al, polars, polarsA] = interpolate_alpha(Alpha, AlphaA, polars, polarsA, Alpha0);
    polars_dA = (polarsA - polars)/(3);
    CldA0 = polars_dA(4);

    [polars, Alpha, alpha, txt] = read_data_static(filenamea.nominal);
    [polarsAS, AlphaAS, alphaAS, txtAS] = read_data_static(filenamea.AS);
    [Al, al, polars, polarsAS] = interpolate_alpha(Alpha, AlphaAS, polars, polarsAS, Alpha0);
    polars_dAS = (polarsAS - polars)/(3);
    CLdAS0 = polars_dAS(1);
    CmdAS0 = polars_dAS(5);
    CDdAS0 = polars_dAS(2);

    [polars, Alpha, alpha, txt] = read_data_static(filenamea.nominal);
    [polarsES, AlphaES, alphaES, txtES] = read_data_static(filenamea.ES);
    [Al, al, polars, polarsES] = interpolate_alpha(Alpha, AlphaES, polars, polarsES, Alpha0);
    polars_dES = (polarsES - polars)/(3);
    CldES0 = polars_dES(4);
    
    [CLdE0, CDdE0, CmdE0; 
     CYdR0, CndR0, CldR0; 
     CldA0, 0, 0; 
     CLdAS0, CDdAS0, CmdAS0; 
     CldES0, 0, 0]
    % Tests
%     chord = 0.35;
%     tail_arm = 1.070;
%     [CLdE0, CmdE0*chord/tail_arm]
%     
%     span = 2.910;
%     tail_span = 0.570;
%     [CLdE0, CldES0*span/(tail_span/2)]

    % extend to 90deg
    CLq = @(A) CLq0 + 0*A;
    Cmq = @(A) Cmq0 + 0*A;
    CYb = @(A) cos(deg2rad(A))*(CYb0/cos(deg2rad(Alpha0)));
    CYr = @(A) CYr0 + 0*A;
    CYp = @(A) CYp0 + 0*A;
    Clp = @(A) Clp0 + 0*A;
    Clb = @(A) cos(deg2rad(A))*(Clb0/cos(deg2rad(Alpha0)));
    Clr = @(A) Clr0 + 0*A;
    Cnb = @(A) cos(deg2rad(A))*(Cnb0/cos(deg2rad(Alpha0)));
    Cnr = @(A) Cnr0 + 0*A;
    Cnp = @(A) Cnp0 + 0*A;

    CLdE = @(A) cos(deg2rad(A))*(CLdE0/cos(deg2rad(Alpha0)));
    CmdE = @(A) cos(deg2rad(A))*(CmdE0/cos(deg2rad(Alpha0)));
    CDdE = @(A) cos(deg2rad(A))*(CDdE0/cos(deg2rad(Alpha0)));
    CYdR = @(A) cos(deg2rad(A))*(CYdR0/cos(deg2rad(Alpha0)));
    CldR = @(A) cos(deg2rad(A))*(CldR0/cos(deg2rad(Alpha0)));
    CndR = @(A) cos(deg2rad(A))*(CndR0/cos(deg2rad(Alpha0)));
    CldA = @(A) cos(deg2rad(A))*(CldA0/cos(deg2rad(Alpha0)));
    CLdAS = @(A) cos(deg2rad(A))*(CLdAS0/cos(deg2rad(Alpha0)));
    CDdAS = @(A) cos(deg2rad(A))*(CDdAS0/cos(deg2rad(Alpha0)));
    CmdAS = @(A) cos(deg2rad(A))*(CmdAS0/cos(deg2rad(Alpha0)));
    CldES = @(A) cos(deg2rad(A))*(CldES0/cos(deg2rad(Alpha0)));

    function [Al, al, polars, polarsE] = interpolate_alpha(Alpha, AlphaE, polars, polarsE, Alpha0)
        Al = Alpha0;
        polars = interp1(Alpha, polars, Alpha0);
        polarsE = interp1(AlphaE, polarsE, Alpha0);
        al = deg2rad(Al);
    end

    function [polars, Alpha, alpha, txt] = read_data_static(filename)
        delimiterIn = ' ';
        headerlinesIn = 8;
        A = importdata(filename,delimiterIn,headerlinesIn);
        polars = A.data;

        Alpha = polars(:,1);
        alpha = deg2rad(Alpha);
        polars = polars(:, [3,6,7,8,9,10]);
        txt = A.textdata;
        txt = txt(end,[3,6,7,8,9,10]);
    end

end


function [] = load_propeller_thrust_model(Rp1, Rp2, Rp3, Rp4, ep1, ep2, ep3, ep4, np1_max, np2_max, np3_max, np4_max, filenamep)
    %% Defines the Thrust model and the Slipstream velocity model
    
    global T1 T2 T3 T4 Q1 Q2 Q3 Q4 Vp1 Vp2 Vp3 Vp4 dP1min dP2min dP3min dP4min P1 P2 P3 P4
    [R_inv] = Rotation_and_Euler_Matrices(ep1);
    [T, Q, Vp, dP1min, P1] = fit_thrust_fn(filenamep.p1, Rp1*2, np1_max, 1);
    T1 = @(V,n) R_inv*[T(V,n); 0; 0];
    Q1 = @(V,n) R_inv*[-Q(V,n);0;0]; % reversed direction
    Vp1 = @(V,n) R_inv*[Vp(V,n); 0; 0];
    [R_inv] = Rotation_and_Euler_Matrices(ep2);
    [T, Q, Vp, dP2min, P2] = fit_thrust_fn(filenamep.p2, Rp2*2, np2_max, 0);
    T2 = @(V,n) R_inv*[T(V,n); 0; 0];
    Q2 = @(V,n) R_inv*[Q(V,n);0;0];
    Vp2 = @(V,n) R_inv*[Vp(V,n); 0; 0];
    [R_inv] = Rotation_and_Euler_Matrices(ep3);
    [T, Q, Vp, dP3min, P3] = fit_thrust_fn(filenamep.p3, Rp3*2, np3_max, 1);
    T3 = @(V,n) R_inv*[T(V,n); 0; 0];
    Q3 = @(V,n) R_inv*[Q(V,n);0;0];
    Vp3 = @(V,n) R_inv*[Vp(V,n); 0; 0];
    [R_inv] = Rotation_and_Euler_Matrices(ep4);
    [T, Q, Vp, dP4min, P4] = fit_thrust_fn(filenamep.p4, Rp4*2, np4_max, 0);
    T4 = @(V,n) R_inv*[T(V,n); 0; 0];
    Q4 = @(V,n) R_inv*[Q(V,n);0;0];
    Vp4 = @(V,n) R_inv*[Vp(V,n); 0; 0];
end

function [T, Q, Vp, dPmin, P, ETA_fit] = fit_thrust_fn(filep, D, npmax, plots)
    global rho fdec

    A = importdata(filep);
    polars = A.data;
    J = polars(:,1);
    CT = polars(:,2);
    CP = polars(:,3);
    ETA = polars(:,4);
    X = fliplr(vander(J));
    X = X(:, 1:4);

    %% CT thrust coef
    y = CT;
    a = (X'*X)\(X'*y);
    % cT = @(V,n) max((a(3)*(V./(n*D)).^2 + a(2)*(V./(n*D)) + a(1)), min(CT));
    % cT = @(V,n) max((a(3)*(V./(n*D)).^2 + a(2)*(V./(n*D)) + a(1)), 0);
    cT = @(V,n) max((a(4)*(V./(n*D)).^3 + a(3)*(V./(n*D)).^2 + a(2)*(V./(n*D)) + a(1)), 0).*((V./(n*D))<2); % cubic fit, trimmed to min 0 and checked against n->0
    %cT = @(V,n) max((a(3)*(V./(n*D)).^2 + a(2)*(V./(n*D)) + a(1)), 0);
    T = @(V,n) cT(V,n)*rho.*(n.^2).*(D.^4);

    %% CP power coef
    y = CP;
    a = (X'*X)\(X'*y);
    cP = @(V,n) max((a(4)*(V./(n*D)).^3 + a(3)*(V./(n*D)).^2 + a(2)*(V./(n*D)) + a(1)), 0).*((V./(n*D))<2);
    %cP = @(V,n) ((a(3)*(V./(n*D)).^2 + a(2)*(V./(n*D)) + a(1)))*((n/npmax)>dPmin(V));
    P = @(V,n) cP(V,n)*rho.*(n.^3)*D^5;

    %% CQ torque coef
    CQ = CP/(2*pi);
    y = CQ;
    a = (X'*X)\(X'*y);
    cQ = @(V,n) max((a(4)*(V./(n*D)).^3 + a(3)*(V./(n*D)).^2 + a(2)*(V./(n*D)) + a(1)), 0).*((V./(n*D))<2);
    Q = @(V,n) cQ(V,n)*rho.*(n.^2)*D^5;
    
    %% ETA 
    ETA_fit = @(V,n) max(interp1(J, ETA, V./(n*D), 'spline', 'extrap'), 0).*((V./(n*D))<2);

    %% Slipstream velo
    A = (pi*D^2)/4;
    vp = @(V,n) sqrt(2*T(V,n)/(rho*A) + V.^2);
    Vp = @(V,n) fdec(V)*(vp(V,n) - V);

    %% Minimum rpm to have some nonzero thrust
    dPmin = @(V) fzero(@(n) T(V,n) - 0.001, [0,2*npmax])/npmax;

    %% Plots
    if plots
        if contains(filep, 'wing') 
            stri = sprintf('wing prop');
        elseif contains(filep, 'front')
            stri = sprintf('front prop');
        else
            stri = sprintf('extrap');
        end
        n = 5000/60;
        Vs = 0:0.5:40;
        figure(109)
        plot(J, CT, 'r*', 'HandleVisibility', 'Off');
        hold on;
        plot(Vs./(n*D), cT(Vs,n), 'DisplayName', stri);
        grid minor;
        legend();
        xlabel('J [-]');
        ylabel('CT [-]');

        figure(110)
        plot(J, CQ, 'r*', 'HandleVisibility', 'Off');
        hold on;
        plot(Vs./(n*D), cQ(Vs,n), 'DisplayName', stri);
        grid minor;
        legend();
        xlabel('J [-]');
        ylabel('CQ [-]');

        figure(111)
        plot(J, CP, 'r*', 'HandleVisibility', 'Off');
        hold on;
        plot(Vs./(n*D), cP(Vs,n), 'DisplayName', stri);
        grid minor;
        legend();
        xlabel('J [-]');
        ylabel('CP [-]');

        figure(112)
        plot(J, ETA, 'r*', 'HandleVisibility', 'Off');
        hold on;
        plot(Vs./(n*D), ETA_fit(Vs,n), 'DisplayName', stri);
        grid minor;
        legend();
        xlabel('J [-]');
        ylabel('ETA [-]');
        
        figure(113)
        hold on;
        ns = [0:1000:10000]/60;
        for i=1:length(ns)
            Ts = zeros(length(Vs),1);
            for j = 1:length(Vs)
                Ts(j) = vp(Vs(j),ns(i));
            end
            plot(Vs, Ts, '--');
            text(Vs(1), Ts(1), sprintf(' n=%2.0f',ns(i)*60));
        end
        grid minor;
        xlabel('V [m/s]');
        ylabel('Vss [m/s]');
        title('Slipstream velocity(V): n from 0 to 10000 rpm');
    end
end

function [] = load_Custer_Channel_lift_effect(filenamecc)
    global fcc
    A = importdata(filenamecc.lift); % from Wind Tunnel data
    A = A.data;
    Vs = A(:,1);
    fs = A(:,2);
    
    fcc = @(V) interp1(Vs, fs, V, 'linear', 'extrap');
end

function [] = load_slipstream_decay_factor()
    global fdec Vcruise
    Vs = [0, Vcruise, 2*Vcruise];
    fs = [0.8, 1.8, 1.8];
    fdec = @(V) interp1(Vs, fs, V, 'linear', 'extrap');
end


function [] = plot_aero_derivatives()
    global CLq Cmq CYb CYr Clp Cnb Cnr
    global CLdE CDdE CmdE CYdR CldR CndR CldA CLdAS CDdAS CmdAS CldES
    Alpha = linspace(0, 90, 500);
    
    figure(201)
    subplot(3,3,3);
    plot(Alpha, CLq(Alpha));
    grid minor;
    title('CLq')
    subplot(3,3,6);
    plot(Alpha, 0*Cmq(Alpha));
    grid minor;
    title('CDq')
    subplot(3,3,9);
    plot(Alpha, Cmq(Alpha));
    grid minor;
    title('Cmq')
    
    figure(202)
    subplot(3,3,1);
    plot(Alpha, CYb(Alpha));
    grid minor;
    title('CYb')
    subplot(3,3,3)
    plot(Alpha, CYr(Alpha));
    grid minor;
    title('CYr')
    subplot(3,3,5);
    plot(Alpha, Clp(Alpha));
    grid minor;
    title('Clp')
    subplot(3,3,7);
    plot(Alpha, Cnb(Alpha));
    grid minor;
    title('Cnb')
    subplot(3,3,9);
    plot(Alpha, Cnr(Alpha));
    grid minor;
    title('Cnr')
    
    figure(204)
    subplot(6,4,3)
    plot(Alpha, CLq(Alpha));
    grid minor;
    title('CLq')
    subplot(6,4,11)
    plot(Alpha, 0*Cmq(Alpha));
    grid minor;
    title('CDq')
    subplot(6,4,19)
    plot(Alpha, Cmq(Alpha));
    grid minor;
    title('Cmq');
    subplot(6,4,5)
    plot(Alpha, CYb(Alpha));
    grid minor;
    title('CYb')
    subplot(6,4,8)
    plot(Alpha, CYr(Alpha));
    grid minor;
    title('CYr')
    subplot(6,4,14)
    plot(Alpha, Clp(Alpha));
    grid minor;
    title('Clp')
    subplot(6,4,21)
    plot(Alpha, Cnb(Alpha));
    grid minor;
    title('Cnb')
    subplot(6,4,24)
    plot(Alpha, Cnr(Alpha));
    grid minor;
    title('Cnr')
    
%     subplot(6,4,13)
%     plot(Alpha, 0*CLq(Alpha));
%     grid minor;
%     title('Clb')
%     subplot(6,4,6)
%     plot(Alpha, 0*CLq(Alpha));
%     grid minor;
%     title('CYp')
%     subplot(6,4,22)
%     plot(Alpha, 0*CLq(Alpha));
%     grid minor;
%     title('Cnp')
%     subplot(6,4,16)
%     plot(Alpha, 0*CLq(Alpha));
%     grid minor;
%     title('Clr')
    
    
    
    
    figure(203)
    subplot(6,5,1);
    plot(Alpha, CLdE(Alpha));
    grid minor;
    title('   CL_{uE}')
    subplot(6,5,4)
    plot(Alpha, CLdAS(Alpha));
    grid minor;
    title('   CL_{uAS}');
    subplot(6,5,21);
    plot(Alpha, CmdE(Alpha));
    grid minor;
    title('   Cm_{uE}');
    subplot(6,5,24);
    plot(Alpha, CmdAS(Alpha));
    grid minor;
    title('   Cm_{uAS}')
    subplot(6,5,7);
    plot(Alpha, CYdR(Alpha));
    grid minor;
    title('   CY_{uR}');
    subplot(6,5,17);
    plot(Alpha, CldR(Alpha));
    grid minor;
    title('   Cl_{uR}');
    subplot(6,5,18);
    plot(Alpha, CldA(Alpha));
    grid minor;
    title('   Cl_{uA}');
    subplot(6,5,20);
    plot(Alpha, CldES(Alpha));
    grid minor;
    title('   Cl_{uES}');
    subplot(6,5,27);
    plot(Alpha, CndR(Alpha));
    grid minor;
    title('   Cn_{uR}');
    subplot(6,5,11);
    plot(Alpha, CDdE(Alpha));
    grid minor;
    title('   CD_{uE}');
    subplot(6,5,14);
    plot(Alpha, CDdAS(Alpha));
    grid minor;
    title('   CD_{uAS}')
    
%     subplot(6,5,11);
%     plot(Alpha, 0*CLdE(Alpha));
%     grid minor;
%     title(' CD_{uE}')
%     subplot(6,5,28);
%     plot(Alpha, 0*CLdE(Alpha));
%     grid minor;
%     title(' Cn_{uA}')
%     subplot(6,5,14);
%     plot(Alpha, 0*CLdE(Alpha));
%     grid minor;
%     title(' CD_{uAS}')
%     subplot(6,5,30);
%     plot(Alpha, 0*CLdE(Alpha));
%     grid minor;
%     title(' Cl_{uES}')
    
end

function [] = plot_thrust_model()
    global T4 Q4 Vp4 dP4min P4
    
    dir_ind = 1; % index of the vector direcion in which this happens
    
    V = linspace(0, 40, 100);
    n = linspace(0, 10000, 100)/60;
    
    figure(100)
    plot(V, V*0, 'k--');
    hold on;
    ns = [0:1000:10000]/60;
    for i=1:length(ns)
        Ts = zeros(length(V),1);
        for j = 1:length(V)
            T = T4(V(j),ns(i));
            Ts(j) = T(dir_ind);
        end
        plot(V, Ts, '--');
        text(V(1), Ts(1), sprintf(' n=%2.0f',ns(i)*60));
    end
    grid minor;
    xlabel('V [m/s]');
    ylabel('T [N]');
    title('Thrust(V): n from 0 to 10000 rpm');

    figure(101)
    hold on;
    Vs = 0:5:40;
    for i=1:length(Vs)
        Ts = zeros(length(V),1);
        for j = 1:length(n)
            T = T4(Vs(i),n(j));
            Ts(j) = T(dir_ind);
        end
        plot(n*60, Ts, '--');
        text(n(end)*60, Ts(end), sprintf('V=%2.0f',Vs(i)));
    end
    grid minor;
    xlabel('n [rpm]');
    ylabel('T [N]');
    title('Thrust(n): V from 0 to 40 m/s');
    
    figure(102)
    hold on;
    ns = [0:1000:10000]/60;
    for i=1:length(ns)
        Ts = zeros(length(V),1);
        for j = 1:length(V)
            T = Vp4(V(j),ns(i));
            Ts(j) = T(dir_ind);
        end
        plot(V, Ts, '--');
        text(V(1), Ts(1), sprintf(' n=%2.0f',ns(i)*60));
    end
    grid minor;
    xlabel('V [m/s]');
    ylabel('Vss [m/s]');
    title('Slipstream velocity(V): n from 0 to 10000 rpm');

    figure(103)
    hold on;
    Vs = 0:5:40;
    for i=1:length(Vs)
        Ts = zeros(length(V),1);
        for j = 1:length(n)
            T = Vp4(Vs(i),n(j));
            Ts(j) = T(dir_ind);
        end
        plot(n*60, Ts, '--');
        text(n(end)*60, Ts(end), sprintf('V=%2.0f',Vs(i)));
    end
    grid minor;
    xlabel('n [rpm]');
    ylabel('Vss [m/s]');
    title('Slipstream velocity(n): V from 0 to 40 m/s');
    
    figure(104)
    hold on;
    for i = 1:length(V)
        plot(V(i), dP4min(V(i)), 'b.', 'MarkerSize', 10);
    end
    xlabel('V');
    ylabel('dPmin');
    title('Minimal dP=np/npmax to produce thrust');
    grid minor;
    
    
    figure(105)
    plot(V, V*0, 'k--');
    hold on;
    ns = [0:1000:10000]/60;
    for i=1:length(ns)
        Ts = zeros(length(V),1);
        for j = 1:length(V)
            T = P4(V(j),ns(i));
            Ts(j) = T;
        end
        plot(V, Ts, '--');
        text(V(1), Ts(1), sprintf(' n=%2.0f',ns(i)*60));
    end
    grid minor;
    xlabel('V [m/s]');
    ylabel('P [W]');
    title('Power(V): n from 0 to 10000 rpm');

    figure(106)
    hold on;
    Vs = 0:5:40;
    for i=1:length(Vs)
        Ts = zeros(length(V),1);
        for j = 1:length(n)
            T = P4(Vs(i),n(j));
            Ts(j) = T;
        end
        plot(n*60, Ts, '--');
        text(n(end)*60, Ts(end), sprintf('V=%2.0f',Vs(i)));
    end
    grid minor;
    xlabel('n [rpm]');
    ylabel('P [W]');
    title('Power(n): V from 0 to 40 m/s');
    
    figure(107)
    plot(V, V*0, 'k--');
    hold on;
    ns = [0:1000:10000]/60;
    for i=1:length(ns)
        Ts = zeros(length(V),1);
        for j = 1:length(V)
            T = Q4(V(j),ns(i));
            Ts(j) = T(dir_ind);
        end
        plot(V, Ts, '--');
        text(V(1), Ts(1), sprintf(' n=%2.0f',ns(i)*60));
    end
    grid minor;
    xlabel('V [m/s]');
    ylabel('Q [N.m]');
    title('Torque(V): n from 0 to 10000 rpm');

    figure(108)
    hold on;
    Vs = 0:5:40;
    for i=1:length(Vs)
        Ts = zeros(length(V),1);
        for j = 1:length(n)
            T = Q4(Vs(i),n(j));
            Ts(j) = T(dir_ind);
        end
        plot(n*60, Ts, '--');
        text(n(end)*60, Ts(end), sprintf('V=%2.0f',Vs(i)));
    end
    grid minor;
    xlabel('n [rpm]');
    ylabel('Q [N.m]');
    title('Torque(n): V from 0 to 40 m/s');

end

function [] = plot_decay_factor()
    global fdec
    figure(400)
    Vs = linspace(0, 50, 100);
    plot(Vs, fdec(Vs), 'b-');
    grid minor;
    title('Slipstream induced velocity decay factor');
    xlabel('V [m/s]');
    ylabel('fdec');
    ylim([0,2]);
end

function [] = plot_Custer_Channel_lift_effect(filenamecc)
    
    global fcc
    figure(300)
    Vs = linspace(0, 15, 40);
    plot(Vs, fcc(Vs), 'b-', 'DisplayName', 'interpolated');
    grid minor;
    title('Custer channel upward-thrust effect: \DeltaZ = -fcc.\DeltaX');
    xlabel('V [m/s]');
    ylabel('fcc');
    
    A = importdata(filenamecc.lift);
    polars = A.data;
    Vs = polars(:,1);
    fs = polars(:,2);
    hold on;
    plot(Vs, fs, 'r*', 'MarkerSize', 15, 'DisplayName', 'measured');
    legend('location', 'best');
    ylim([0,1]);
end
