function Load_params(cs)
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
            % rho = 1.179; % to match X-Plane
            g = 9.81; % [m/s^2]
            
            %% Aerodynamics
            Vcruise = 25; % cruise speed [m/s]
            % static and stability
            XCP = -0.35/4;
            filenamestatic.CL = fullfile('Aero_data_files4', 'CL.txt');
            filenamestatic.CD = fullfile('Aero_data_files4', 'CD.txt');
            filenamea.SD = fullfile('Aero_data_files4', 'T7-raw-out0.txt');
            % control
            udA_max = 40; % 39
            udE_max = 30; % [deg] max deflection of control surfaces %27
            udR_max = 40; % 17.5
            udAS_max = 10;
            udES_max = 10; % 9.5
            x_EL = [-1.07, -0.26, 0.00];
            x_RL = [-1.07, -0.57, -0.20];
            x_AL = [-0.20, -0.99, 0.00];
            filenamea.E0 = fullfile('Aero_data_files4', 'T1-20_0 m_s-VLM2-16_0kg-x39_1mm-dE00.txt');
            filenamea.R0 = fullfile('Aero_data_files4', 'T1-20_0 m_s-VLM2-16_0kg-x39_1mm-dR00.txt');
            filenamea.A0 = fullfile('Aero_data_files4', 'T1-20_0 m_s-VLM2-16_0kg-x39_1mm-dA00.txt');
            filenamea.AS0 = fullfile('Aero_data_files4', 'T1-20_0 m_s-VLM2-16_0kg-x39_1mm-dAS00.txt');
            filenamea.E3 = fullfile('Aero_data_files4', 'T1-20_0 m_s-VLM2-16_0kg-x39_1mm-dE03.txt');
            filenamea.R3 = fullfile('Aero_data_files4', 'T1-20_0 m_s-VLM2-16_0kg-x39_1mm-dR03.txt');
            filenamea.A3 = fullfile('Aero_data_files4', 'T1-20_0 m_s-VLM2-16_0kg-x39_1mm-dA03.txt');
            filenamea.AS3 = fullfile('Aero_data_files4', 'T1-20_0 m_s-VLM2-16_0kg-x39_1mm-dAS03.txt');
            
            %% Custer effect
            filenamecc.lift = fullfile('Aero_data_files4', 'custer_effect.txt');
            filenamecc.stall.CL = fullfile('Aero_data_files4', 'CLCuster.txt');
            filenamecc.stall.CD = fullfile('Aero_data_files4', 'CDCuster.txt');
            
            %% Propellers
            filenamep.p1 = fullfile('Aero_data_files4', 'my_prop_wing.txt');
            filenamep.p2 = fullfile('Aero_data_files4', 'my_prop_wing.txt');
            filenamep.p3 = fullfile('Aero_data_files4', 'my_prop_front.txt');
            filenamep.p4 = fullfile('Aero_data_files4', 'my_prop_front.txt'); % just for plotting
            
            mp1 = 0.072; % [kg] % left
            Rp1 = 0.356/2; % [m] % 0.34 the shorter one but in WT it made no difference
            ep1 = [0,deg2rad(8),0]; % [rad]
            xp1 = [0*-0.35,-0.33,0]; % [m]
            np1_max = 9000/60; % [rps]
            mh1 = 0.15; % [kg] elecetric motor rotor mass
            Rh1 = 0.058/2; % [m] electric motor rotor radius
            
            mp2 = 0.072; % right
            Rp2 = 0.356/2;
            ep2 = [0,deg2rad(8),0];
            xp2 = [0*-0.35,0.33,0];
            np2_max = 9000/60;
            mh2 = 0.15; % [kg] elecetric motor rotor mass
            Rh2 = 0.058/2; % [m] electric motor rotor radius
            
            mp3 = 0.09; % front
            Rp3 = 0.457/2;
            ep3 = [0,0,0];
            xp3 = [0.8,0,0];
            np3_max = 8000/60;
            mh3 = 0.15; % [kg] front engine hub mass
            Rh3 = 0.058/2; % [m] front engine hub radius
            
            mp4 = 0; % for plotting
            Rp4 = 0.17;
            ep4 = [0,deg2rad(0),0];
            xp4 = [0,0,0];
            np4_max = 9000/60;
            mh4 = 0.15; 
            Rh4 = 0.058/2; 
        
        otherwise
            fprintf('Error: Case %d not defined.\n', cs);
    end
    
    
    %%% Processing
    % Read the aerodynamic derivatives data
    load_aerodynamic_derivatives(filenamea, x_EL, x_RL, x_AL);
    
    % Read the static aerodynamic coefs with the custer stall effects
    Aero_static_fits(filenamestatic, filenamecc);
    
    
    minv = 1/m;
    Iinv = inv(I);
    
    % Propeller inertial terms
    [Ip1] = Mom_Inertia_estim_prop(mp1, Rp1, ep1, mh1, Rh1);
    [Ip2] = Mom_Inertia_estim_prop(mp2, Rp2, ep2, mh2, Rh2);
    [Ip3] = Mom_Inertia_estim_prop(mp3, Rp3, ep3, mh3, Rh3);
    [Ip4] = Mom_Inertia_estim_prop(mp4, Rp4, ep4, mh4, Rh4);
    
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

function [Ip] = Mom_Inertia_estim_prop(m, R, ep, mh, Rh)
    % https://www.mh-aerotools.de/airfoils/prop_precession_english.htm
    % Three rods of a mass M=m/3 spinning around their end. 
    % One rod Ix = 1/3*ML^2
    % Three rods of m Ix = 1/3mL^2
    % Spinning disc is better as an upper bound Ix = 1/2mR^2
    
    % Also hub I added, for example the electric motor rotor is quite heavy
    % Spinning disc Ix = 1/2mR^2
    Ix = 0.5*m*R^2 + 0.5*mh*Rh^2;
    Ip = zeros(3,3);
    Ip(1,1) = Ix;
    [~,~,~,~, R, ~,~,~] = Rotation_and_Euler_Matrices(ep);
    % Ip = R'*Ip*R; % prop to the body frame
    % Ip = Ip*R'; % so that in Hp = Ip.omegap = R'IpR . R'omegap = R'IpRR' . omegap you can use the omegap in the prop frame so [omgx; 0; 0] vector
    Ip = R'*Ip; % result of both up together. Not really Ip in the body frame but takes care of the transformation when multiplied by omega
end


function [] = load_aerodynamic_derivatives(filenamea, x_EL, x_RL, x_AL)
    
    global CLq Cmq CYb CYr Clp Cnb Cnr
    global CYp Clb Clr Cnp
    global CLdEL CDdEL CmdEL CldEL CYdRL CldRL CndRL CldAL CLdAL CDdAL CmdAL
    global c b
        
    %% Aerodynamic Derivatives 
    delimiterIn = ',';
    headerlinesIn = 1;
    A = importdata(filenamea.SD,delimiterIn,headerlinesIn);
    polars = A.data;
    Alpha0 = 0;
    polars = interp1(polars(:,1), polars, Alpha0);
    CLq0 = polars(31); 
    Cmq0 = polars(34); 
    CYb0 = -polars(35); 
    CYp0 = 1*polars(36); % very small
    CYr0 = polars(37); 
    Clb0 = -1*polars(38); % small
    Clp0 = polars(39); 
    Clr0 = 1*polars(40); % small
    Cnb0 = -polars(41); 
    Cnp0 = 1*polars(42); % adverse yaw, very small
    Cnr0 = polars(43); 
    
    if not(isempty(strfind(filenamea.SD, 'Aero_data_files4')))
        Cnr0 = 1*Cnr0; % you can modify a coefs value here for case2 comparison
    end

    %% Control Derivatives
    [polars, Alpha, alpha, txt] = read_data_static(filenamea.E0);
    [polarsE, AlphaE, alphaE, txtE] = read_data_static(filenamea.E3);
    [Al, al, polars, polarsE] = interpolate_alpha(Alpha, AlphaE, polars, polarsE, Alpha0);
    polars_dE = (polarsE - polars)/(3);
    CLdEL0 = polars_dE(1)/2; % /2 because only left elevator, both used in xflr to avoid doubling the number of elevator tips
    CDdEL0 = polars_dE(2)/2;
    CYdEL0 = 0;
    cM = Cross_Matrix(x_EL)*[-CDdEL0; CYdEL0; -CLdEL0];
    CmdEL0 = cM(2)/c; % non-dim the arms
    CldEL0 = cM(1)/b;

    [polars, Alpha, alpha, txt] = read_data_static(filenamea.R0);
    [polarsR, AlphaR, alphaR, txtR] = read_data_static(filenamea.R3);
    [Al, al, polars, polarsR] = interpolate_alpha(Alpha, AlphaR, polars, polarsR, Alpha0);
    polars_dR = (polarsR - polars)/(3);
    CLdRL0 = 0;
    CDdRL0 = 0;
    CYdRL0 = polars_dR(1); % only left rudder
    cM = Cross_Matrix(x_RL)*[-CDdRL0; CYdRL0; -CLdRL0];
    CldRL0 = cM(1)/b;
    CndRL0 = cM(3)/b; 

    [polars, Alpha, alpha, txt] = read_data_static(filenamea.AS0);
    [polarsAS, AlphaAS, alphaAS, txtAS] = read_data_static(filenamea.AS3);
    [Al, al, polars, polarsAS] = interpolate_alpha(Alpha, AlphaAS, polars, polarsAS, Alpha0);
    polars_dAS = (polarsAS - polars)/(3);
    CLdAL0 = polars_dAS(1)/2; % /2 because left only
    CDdAL0 = polars_dAS(2)/2;
    CYdAL0 = 0;
    cM = Cross_Matrix(x_AL)*[-CDdAL0; CYdAL0; -CLdAL0];
    CldAL0 = cM(1)/b;
    CmdAL0 = cM(2)/c;
    
%     [CLdEL0, CDdEL0, CmdEL0; 
%      CYdRL0, CndRL0, CldRL0; 
%      CldAL0, 0, 0; 
%      CLdAL0, CDdAL0, CmdAL0; 
%      CldEL0, 0, 0]

    % extend to 90deg
    CLq = @(A) CLq0 + 0*A;
    Cmq = @(A) Cmq0 + 0*A;
    CYb = @(A) cosd(A)*(CYb0/cosd(Alpha0));
    CYr = @(A) CYr0 + 0*A;
    CYp = @(A) CYp0 + 0*A;
    Clp = @(A) Clp0 + 0*A;
    Clb = @(A) cosd(A)*(Clb0/cosd(Alpha0));
    Clr = @(A) Clr0 + 0*A;
    Cnb = @(A) cosd(A)*(Cnb0/cosd(Alpha0));
    Cnr = @(A) Cnr0 + 0*A;
    Cnp = @(A) Cnp0 + 0*A;

    CLdEL = @(A) cosd(A)*(CLdEL0/cosd(Alpha0));
    CmdEL = @(A) cosd(A)*(CmdEL0/cosd(Alpha0));
    CDdEL = @(A) cosd(A)*(CDdEL0/cosd(Alpha0));
    CldEL = @(A) cosd(A)*(CldEL0/cosd(Alpha0));
    CYdRL = @(A) cosd(A)*(CYdRL0/cosd(Alpha0));
    CldRL = @(A) cosd(A)*(CldRL0/cosd(Alpha0));
    CndRL = @(A) cosd(A)*(CndRL0/cosd(Alpha0));
    CldAL = @(A) cosd(A)*(CldAL0/cosd(Alpha0));
    CLdAL = @(A) cosd(A)*(CLdAL0/cosd(Alpha0));
    CDdAL = @(A) cosd(A)*(CDdAL0/cosd(Alpha0));
    CmdAL = @(A) cosd(A)*(CmdAL0/cosd(Alpha0));
    

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
    Q1 = @(V,n) R_inv*[-Q(V,n);0;0]; % -Q to reverse direction if counter-spinning
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
        hold on;
        %plot(J, CT, 'r*', 'HandleVisibility', 'Off');
        %hold on;
        plot(Vs./(n*D), cT(Vs,n), 'DisplayName', stri);
        grid minor;
        legend();
        xlabel('$\frac{V_p}{n_pD_p}$ [-]', 'Interpreter', 'latex', 'FontSize', 15);
        ylabel('$C_T^{\,p}$ [-]', 'Interpreter', 'latex');

        figure(110)
        hold on;
        %plot(J, CQ, 'r*', 'HandleVisibility', 'Off');
        %hold on;
        plot(Vs./(n*D), cQ(Vs,n), 'DisplayName', stri);
        grid minor;
        legend();
        xlabel('$\frac{V_p}{n_pD_p}$ [-]', 'Interpreter', 'latex', 'FontSize', 15);
        ylabel('$C_Q^{\,p}$ [-]', 'Interpreter', 'latex');

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
        
%         figure(114)
%         hold on;
%         plot(Vs./(n*D), ETA_fit(Vs,n), 'DisplayName', stri);
%         plot(Vs./(n*D), (Vs./(n*D)).*cT(Vs,n)./cP(Vs,n), 'DisplayName', stri);
%         grid minor;
%         legend();
%         xlabel('J [-]');
%         ylabel('ETA [-]');
%         title('check for a formula');
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
    global CYp Clb Clr Cnp
    global CLdEL CDdEL CmdEL CldEL CYdRL CldRL CndRL CldAL CLdAL CDdAL CmdAL
    Alpha = linspace(0, 91, 500);
    
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
    
    ticks = [0,45,90];
    tickslbl = {0,'- \alpha -',90};
    figure(204)
    subplot(6,4,3)
    plot(Alpha, CLq(Alpha));
    grid minor;
    title('CLq')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(CLq(0),1)]);
%     subplot(6,4,11)
%     plot(Alpha, 0*Cmq(Alpha));
%     grid minor;
%     title('CDq')
    subplot(6,4,19)
    plot(Alpha, Cmq(Alpha));
    grid minor;
    title('Cmq');
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(Cmq(0),1)]);
    subplot(6,4,5)
    plot(Alpha, CYb(Alpha));
    grid minor;
    title('CY\beta')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([0, round(CYb(0), 2)]);
    subplot(6,4,8)
    plot(Alpha, CYr(Alpha));
    grid minor;
    title('CYr')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(CYr(0), 2)]);
    subplot(6,4,6)
    plot(Alpha, CYp(Alpha));
    grid minor;
    title('CYp')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(CYp(0), 2)]);
    subplot(6,4,13)
    plot(Alpha, Clb(Alpha));
    grid minor;
    title('Cl\beta')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([0, round(Clb(0), 2)]);
    subplot(6,4,14)
    plot(Alpha, Clp(Alpha));
    grid minor;
    title('Clp')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(Clp(0), 2)]);
    subplot(6,4,16)
    plot(Alpha, Clr(Alpha));
    grid minor;
    title('Clr')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(Clr(0), 2)]);
    subplot(6,4,21)
    plot(Alpha, Cnb(Alpha));
    grid minor;
    title('Cn\beta')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(Cnb(0), 2), 0]);
    subplot(6,4,22)
    plot(Alpha, Cnp(Alpha));
    grid minor;
    title('Cnp')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(Cnp(0), 2)]);
    subplot(6,4,24)
    plot(Alpha, Cnr(Alpha));
    grid minor;
    title('Cnr')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(Cnr(0), 2)]);
    
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
    subplot(6,3,1);
    plot(Alpha, CLdAL(Alpha));
    grid minor;
    title('   CL_{uAL}')
    subplot(6,3,7);
    plot(Alpha, CDdAL(Alpha));
    grid minor;
    title('   CD_{uAL}')
    subplot(6,3,13);
    plot(Alpha, CmdAL(Alpha));
    grid minor;
    title('   Cm_{uAL}')
    subplot(6,3,10);
    plot(Alpha, CldAL(Alpha));
    grid minor;
    title('   Cl_{uAL}')
    
    subplot(6,3,2);
    plot(Alpha, CLdEL(Alpha));
    grid minor;
    title('   CL_{uEL}')
    subplot(6,3,8);
    plot(Alpha, CDdEL(Alpha));
    grid minor;
    title('   CD_{uEL}')
    subplot(6,3,14);
    plot(Alpha, CmdEL(Alpha));
    grid minor;
    title('   Cm_{uEL}')
    subplot(6,3,11);
    plot(Alpha, CldEL(Alpha));
    grid minor;
    title('   Cl_{uEL}')
    
    subplot(6,3,6);
    plot(Alpha, CYdRL(Alpha));
    grid minor;
    title('   CY_{uRL}');
    subplot(6,3,12);
    plot(Alpha, CldRL(Alpha));
    grid minor;
    title('   Cl_{uRL}');
    subplot(6,3,18);
    plot(Alpha, CndRL(Alpha));
    grid minor;
    title('   Cn_{uRL}');  
    
    
    
    
    
    
    
    ticks = [0,45,90];
    tickslbl = {0,'- \alpha -',90};
    figure(205)
    set(gcf, 'Position', [961   152   120   588]);
    subplot(6,1,1);
    plot(Alpha, CLdAL(Alpha));
    grid minor;
    title('CL_{\deltaAP}', 'horizontalAlignment', 'left')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([0, round(CLdAL(0), 4)]);
    
    subplot(6,1,4);
    plot(Alpha, CldAL(Alpha));
    grid minor;
    title('Cl_{\deltaAP}', 'horizontalAlignment', 'left')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([0, round(CldAL(0), 5)]);
    
    subplot(6,1,5);
    plot(Alpha, CmdAL(Alpha));
    grid minor;
    title('Cm_{\deltaAP}', 'horizontalAlignment', 'left')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(CmdAL(0), 5), 0]);
    
    figure(206)
    set(gcf, 'Position', [961   152   120   588]);
    subplot(6,1,1);
    plot(Alpha, CLdAL(Alpha));
    grid minor;
    title('CL_{\deltaAS}', 'horizontalAlignment', 'left')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([0, round(CLdAL(0), 4)]);
    
    subplot(6,1,4);
    plot(Alpha, -CldAL(Alpha));
    grid minor;
    title('Cl_{\deltaAS}', 'horizontalAlignment', 'left')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(-CldAL(0), 5), 0]);
    
    subplot(6,1,5);
    plot(Alpha, CmdAL(Alpha));
    grid minor;
    title('Cm_{\deltaAS}', 'horizontalAlignment', 'left')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(CmdAL(0), 5), 0]);
    
    
    figure(207)
    set(gcf, 'Position', [961   152   120   588]);
    subplot(6,1,2);
    plot(Alpha, CYdRL(Alpha));
    grid minor;
    title('CY_{\deltaRP}', 'horizontalAlignment', 'left');
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([0, round(CYdRL(0), 4)]);
    
    subplot(6,1,4);
    plot(Alpha, CldRL(Alpha));
    grid minor;
    title('Cl_{\deltaRP}', 'horizontalAlignment', 'left');
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([0, round(CldRL(0), 5)]);
    
    subplot(6,1,6);
    plot(Alpha, CndRL(Alpha));
    grid minor;
    title('Cn_{\deltaRP}', 'horizontalAlignment', 'left');
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(CndRL(0), 5), 0]);
    
    
    figure(208)
    set(gcf, 'Position', [961   152   120   588]);
    subplot(6,1,4);
    plot(Alpha, CldRL(Alpha));
    grid minor;
    title('Cl_{\deltaRS}', 'horizontalAlignment', 'left')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([0, round(CldRL(0), 5)]);
    
    subplot(6,1,2);
    plot(Alpha, CYdRL(Alpha));
    grid minor;
    title('CY_{\deltaRS}', 'horizontalAlignment', 'left')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([0, round(CYdRL(0), 4)]);
    
    subplot(6,1,6);
    plot(Alpha, CndRL(Alpha));
    grid minor;
    title('Cn_{\deltaRS}', 'horizontalAlignment', 'left')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(CndRL(0), 5), 0]);
    
    
    
    figure(209)
    set(gcf, 'Position', [961   152   120   588]);
    ticks = [0,45,90];
    tickslbl = {0,'- \alpha -',90};
    subplot(6,1,1);
    plot(Alpha, CLdEL(Alpha));
    grid minor;
    title('CL_{\deltaEP}', 'horizontalAlignment', 'left');
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([0, round(CLdEL(0), 4)]);
    
    subplot(6,1,3);
    plot(Alpha, CDdEL(Alpha));
    grid minor;
    title('CD_{\deltaEP}', 'horizontalAlignment', 'left');
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([0, round(CDdEL(0), 5)]);
    
    subplot(6,1,4);
    plot(Alpha, CldEL(Alpha));
    grid minor;
    title('Cl_{\deltaEP}', 'horizontalAlignment', 'left');
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([0, round(CldEL(0), 5)]);
     
    subplot(6,1,5);
    plot(Alpha, CmdEL(Alpha));
    grid minor;
    title('Cm_{\deltaEP}', 'horizontalAlignment', 'left');
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([ceil(CmdEL(0)*1000)/1000, 0]);
    

    
    ticks = [0,45,90];
    tickslbl = {0,'- \alpha -',90};
    figure(210)
    set(gcf, 'Position', [961   152   120   588]);
    subplot(6,1,1);
    plot(Alpha, CLdEL(Alpha));
    grid minor;
    title('CL_{\deltaES}', 'horizontalAlignment', 'left')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([0, round(CLdEL(0), 4)]);
    
    subplot(6,1,3);
    plot(Alpha, CDdEL(Alpha));
    grid minor;
    title('CD_{\deltaES}', 'horizontalAlignment', 'left')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([0, round(CDdEL(0), 5)]);
    
    subplot(6,1,4);
    plot(Alpha, -CldEL(Alpha));
    grid minor;
    title('Cl_{\deltaES}', 'horizontalAlignment', 'left')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([round(-CldEL(0), 5),0]);
    
    subplot(6,1,5);
    plot(Alpha, CmdEL(Alpha));
    grid minor;
    title('Cm_{\deltaES}', 'horizontalAlignment', 'left')
    xlim([0,91]);
    xticks(ticks);
    xticklabels(tickslbl);
    yticks([ceil(CmdEL(0)*1000)/1000, 0]);
    
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
        if i == 1
            plot(n*60, Ts, '--');
        else
            plot(n*60, Ts, '--', 'HandleVisibility', 'off');
        end
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
        if i == 1
            plot(n*60, Ts, '--');
        else
            plot(n*60, Ts, '--', 'HandleVisibility', 'off');
        end
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
        if i == 1
            plot(n*60, Ts, '--');
        else
            plot(n*60, Ts, '--', 'HandleVisibility', 'off');
        end
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
    Vs = linspace(0, 40, 100);
    plot(Vs, fdec(Vs), 'b-');
    grid minor;
    title('Slipstream velocity decay factor');
    xlabel('$V_p\,[\textrm{m/s}]$', 'Interpreter', 'Latex', 'Fontsize', 14);
    ylabel('$f_{dec}^{p}$', 'Interpreter', 'Latex', 'Fontsize', 14);
    ylim([0,2]);
end

function [] = plot_Custer_Channel_lift_effect(filenamecc)
    
    global fcc
    figure(300)
    Vs = linspace(0, 15, 40);
    plot(Vs, fcc(Vs), 'b-', 'DisplayName', 'interpolated');
    grid minor;
    title('Channel-wing upward-thrust factor');
    xlabel('$V_p\,[\textrm{m/s}]$', 'Interpreter', 'Latex', 'Fontsize', 14);
    ylabel('$f_{cw}^{p}\,[-]$', 'Interpreter', 'Latex', 'Fontsize', 14);
    
    A = importdata(filenamecc.lift);
    polars = A.data;
    Vs = polars(:,1);
    fs = polars(:,2);
    hold on;
    plot(Vs, fs, 'r*', 'MarkerSize', 15, 'DisplayName', 'measured');
    legend('location', 'best');
    ylim([0,1]);
end
