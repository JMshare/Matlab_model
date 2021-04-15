function [F, M, hp, Fa, Fg, Ft_components, Fcc, Fcc_stall, Ma, Mt_components, Mcc, Mcc_stall, Fa_tail, Ma_tail, VssL, Vssf] = Forces_and_Moments(y,u)
% Returns the Force vector F of external forces and moment vector M of
% external moments 

    global m S c b 
    global Ip1 Ip2 Ip3 Ip4 xp1 xp2 xp3 xp4
    global g rho
    global CLstatic CDstatic Cmstatic
    global CLq Cmq CYb CYr Clp Cnb Cnr
    global CYp Clb Clr Cnp
    global CLdEL CDdEL CmdEL CldEL CYdRL CldRL CndRL CldAL CLdAL CDdAL CmdAL
    global T1 T2 T3 T4 Q1 Q2 Q3 Q4 P1 P2 P3 P4
    global fcc CL_Custer_stall CD_Custer_stall Cm_Custer_stall
    global U_max
    global Vp1 Vp2 Vp3 Vp4
    
    
    %% Aerodynamic forces
    % Extract states
    p = y(7);
    q = y(8);
    r = y(9);
    epsilon = y(10:12);
    u = U_max*u;
    dAL = -u(1) + u(4);
    dAR =  u(1) + u(4);
    dEL = u(2) - u(5);
    dER = u(2) + u(5);
    dRL = u(3);
    dRR = u(3);
    dP1 = u(6);
    dP2 = u(7);
    dP3 = u(8);
    dP4 = u(9);
    u = inv(U_max)*u;
    
    
    v = y(4:6);
    V1 = v'*Vp1(5,100)/norm(Vp1(5,100)); % corrections of v to axial direction of the prop
    V2 = v'*Vp2(5,100)/norm(Vp2(5,100));
    V3 = v'*Vp3(5,100)/norm(Vp3(5,100));
    V4 = v'*Vp4(5,100)/norm(Vp4(5,100));
    
    % Compute aerodyn states
    [V, Alpha, Beta] = Extract_aero_states(y);
    qbarS = 0.5*rho*(V^2)*S;
    
%     % Compute RPS from the power setting % discontinuity for the trim
%     routine!
%     dP1 = fsolve(@(n)P1(V,n)-1000*u(6), u(6)*U_max(6,6));
%     dP2 = fsolve(@(n)P2(V,n)-1000*u(7), u(6)*U_max(7,7));
%     dP3 = fsolve(@(n)P3(V,n)-1500*u(8), u(6)*U_max(8,8));
%     dP4 = fsolve(@(n)P4(V,n)-1500*u(6), u(6)*U_max(9,9));
    
    % Compute aerodyn states at the tail
    [Alpha_tailL, V_tailL, Alpha_tailR, V_tailR, VssL, VssR, Vssf] = Propeller_slipstream_effect(dP1, dP2, dP3, V1, V2, V3);
    qbarS_tailL = 0.5*rho*(V_tailL^2)*S;
    qbarS_tailR = 0.5*rho*(V_tailR^2)*S;
    
    % Evaluete the aerodynamic functions at the current state and control
    R_V = @(Alpha) Wind_frame_matrix(Alpha);
    C_DF = @(qbars) diag([qbars, qbars, qbars, qbars*b,qbars*c, qbars*b]);
    C_0 = [CLstatic(Alpha); 0; CDstatic(Alpha); 0; Cmstatic(Alpha); 0];
    C_s = [0,          0,          CLq(Alpha), 0;
           CYb(Alpha), CYp(Alpha), 0,          CYr(Alpha);
           0,          0,          0,          0;
           Clb(Alpha), Clp(Alpha), 0,          Clr(Alpha);
           0,          0,          Cmq(Alpha), 0;
           Cnb(Alpha), Cnp(Alpha), 0,          Cnr(Alpha)];
%     C_uW = [CLdAL(Alpha),   CLdAL(Alpha);
%             0,              0;
%             0*CDdAL(Alpha), 0*CDdAL(Alpha); % can't have less drag if negative dA. should be |dA|, but then non-linear. Neglect.
%             CldAL(Alpha),   -CldAL(Alpha);
%             CmdAL(Alpha),   CmdAL(Alpha);
%             0,              0];
%     C_uL =  @(Alpha) ...
%             [CLdEL(Alpha) ,  0;
%              0            ,  CYdRL(Alpha);
%              CDdEL(Alpha) ,  0;
%              CldEL(Alpha),  CldRL(Alpha);
%              CmdEL(Alpha) ,  0;
%              0            ,  CndRL(Alpha)];
%     C_uR =  @(Alpha) ...
%             [CLdEL(Alpha) ,  0;
%              0            ,  CYdRL(Alpha);
%              CDdEL(Alpha) ,  0;
%              -CldEL(Alpha) ,  CldRL(Alpha);
%              CmdEL(Alpha) ,  0;
%              0            ,  CndRL(Alpha)];
%     
%     FM_a = R_V(Alpha)*C_DF(qbarS)*(C_0 + C_s*C_Ds_i(V)*[deg2rad(Beta); p; q; r] + C_uW*[dAL; dAR]);
%     FM_a_tail = ...
%         + 0.5*R_V(Alpha_tailL)*C_DF(qbarS_tailL)*C_uL(Alpha_tailL)*[(Alpha_tailL+dEL); dRL] ... % half in slipstream
%         + 0.5*R_V(Alpha)*C_DF(qbarS)*C_uL(Alpha)*[(Alpha+dEL); dRL] ... % half not
%         + 0.5*R_V(Alpha_tailR)*C_DF(qbarS_tailR)*C_uR(Alpha_tailR)*[(Alpha_tailR+dER); dRR] ...
%         + 0.5*R_V(Alpha)*C_DF(qbarS)*C_uR(Alpha)*[(Alpha+dER); dRR];
    C_uW = [CLdAL(Alpha),   CLdAL(Alpha),   0,              0;
            0,              0,              CYdRL(Alpha),   CYdRL(Alpha);
            0*CDdAL(Alpha), 0*CDdAL(Alpha), 0,              0; % can't have less drag if negative dA. should be |dA|, but then non-linear. Neglect.
            CldAL(Alpha),   -CldAL(Alpha),  CldRL(Alpha),   CldRL(Alpha);
            CmdAL(Alpha),   CmdAL(Alpha),   0,              0;
            0,              0,              CndRL(Alpha),   CndRL(Alpha)];
    C_uL =  @(Alpha) ...
            [CLdEL(Alpha);
             0;
             CDdEL(Alpha);
             CldEL(Alpha);
             CmdEL(Alpha);
             0];
    C_uR =  @(Alpha) ...
            [CLdEL(Alpha);
             0;
             CDdEL(Alpha);
             -CldEL(Alpha);
             CmdEL(Alpha);
             0];
    
    FM_a = R_V(Alpha)*C_DF(qbarS)*(C_0 + C_s*C_Ds_i(V)*[deg2rad(Beta); p; q; r] + C_uW*[dAL; dAR; dRL; dRR]);
    FM_a_tail = ...
        + 0.5*R_V(Alpha_tailL)*C_DF(qbarS_tailL)*C_uL(Alpha_tailL)*[(Alpha_tailL+dEL)] ... % half in slipstream
        + 0.5*R_V(Alpha)*C_DF(qbarS)*C_uL(Alpha)*[(Alpha+dEL)] ... % half not
        + 0.5*R_V(Alpha_tailR)*C_DF(qbarS_tailR)*C_uR(Alpha_tailR)*[(Alpha_tailR+dER)] ...
        + 0.5*R_V(Alpha)*C_DF(qbarS)*C_uR(Alpha)*[(Alpha+dER)];
  
    Fa_tail = FM_a_tail(1:3);
    Ma_tail = FM_a_tail(4:6);
    Fa = FM_a(1:3) + Fa_tail;
    Ma = FM_a(4:6) + Ma_tail;
     
    %% Gravity forces
    [~, ~, ~, ~, R] = Rotation_and_Euler_Matrices(epsilon);
    Fg = R*[0;0;m*g];
    
    %% Thrust forces
    Ft_components = [T1(V1, dP1), T2(V2, dP2), T3(V3, dP3), T4(V4, dP4)];
    Ft = sum(Ft_components, 2);
    Mt_components = [
        Q1(V1, dP1) + Cross_Matrix(xp1)*Ft_components(:,1), ...
        Q2(V2, dP2) + Cross_Matrix(xp2)*Ft_components(:,2), ...
        Q3(V3, dP3) + Cross_Matrix(xp3)*Ft_components(:,3), ...
        Q4(V4, dP4) + Cross_Matrix(xp4)*Ft_components(:,4)];
    Mt = sum(Mt_components, 2);
    
    %% Custer Ducts perpendicular thrust effect
    Tcc1 = Ft_components(:,1)*fcc(V1);
    Fcc1 = [Tcc1(3); Tcc1(2); -Tcc1(1)];
    Tcc2 = Ft_components(:,2)*fcc(V2);
    Fcc2 = [Tcc2(3); Tcc2(2); -Tcc2(1)];
    Fcc = Fcc1 + Fcc2;
    Mcc = Cross_Matrix(xp1)*Fcc1 + Cross_Matrix(xp2)*Fcc2;
    
    %% Custer Ducts stall delay effect
%     pfract = norm(T1(V1, dP1))/norm(T1(V1, U_max(6,6)));
%     if isnan(pfract) % check against T1(V1,Umax)=0
%         pfract = 0;
%     end 
    pfract = u(6);
    CL = 0.5*CL_Custer_stall(V, Alpha)*pfract; % 0.5 because per one wing
    CD = 0.5*CD_Custer_stall(V, Alpha)*pfract;
    Cm = 0.5*Cm_Custer_stall(V, Alpha)*pfract;
    [CZ, CX] = to_body_axes(CL, CD, Alpha);
    Fcc_stall1 = [CX*qbarS; 0; CZ*qbarS];
    Mcc_stall1 = [0; Cm*qbarS*c; 0] + Cross_Matrix(xp1)*Fcc_stall1;
    
%     pfract = norm(T2(V2, dP2))/norm(T2(V2, U_max(7,7)));
%     if isnan(pfract) % check against T1(V2,Umax)=0
%         pfract = 0;
%     end 
    pfract = u(7);
    CL = 0.5*CL_Custer_stall(V, Alpha)*pfract;
    CD = 0.5*CD_Custer_stall(V, Alpha)*pfract;
    Cm = 0.5*Cm_Custer_stall(V, Alpha)*pfract;
    [CZ, CX] = to_body_axes(CL, CD, Alpha);
    Fcc_stall2 = [CX*qbarS; 0; CZ*qbarS];
    Mcc_stall2 = [0; Cm*qbarS*c; 0] + Cross_Matrix(xp2)*Fcc_stall2;
    
    Fcc_stall = Fcc_stall1 + Fcc_stall2;
    Mcc_stall = Mcc_stall1 + Mcc_stall2;
    
    %% Propeller rotational momentum
    omgp1 = dP1*(2*pi); % [rev/s] to [rad/s]
    omgp2 = dP2*(2*pi);
    omgp3 = dP3*(2*pi);
    omgp4 = dP4*(2*pi);
    hp = Ip1*[-omgp1;0;0] + Ip2*[omgp2;0;0] + Ip3*[omgp3;0;0] + 0*Ip4*[omgp4;0;0];
    
    
    
    %% Total
    F = Fa + Fg + Ft + Fcc + Fcc_stall;
    M = Ma + Mt + Mcc + Mcc_stall;
    
    if(any(isnan([F;M])))
        fprintf('Forces evaluation returns NaN!\n');
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [C_Ds_inv] = C_Ds_i(V)
        if V == 0
            C_Ds_inv = diag([1; 0; 0; 0]);
        else
            C_Ds = diag([1; 2*V/b; 2*V/c; 2*V/b]);
            C_Ds_inv = inv(C_Ds);
        end
    end

    function [Alpha_tail1, V_tail1, Alpha_tail2, V_tail2, VssL, VssR, Vssf] = Propeller_slipstream_effect(dP1, dP2, dP3, V1, V2, V3)
        vp3 = Vp3(V3,dP3);
        v_tail0 = v+vp3/4; % /4 cos the area of the tail is 4 times larger than the prop disc area
        Vssf = norm(vp3);
        
        vp1 = 1*Vp1(V1,dP1);
        v_tail1 = (v_tail0+vp1);
        V_tail1 = norm(v_tail1);
        Alpha_tail1 = rad2deg(atan2(v_tail1(3), v_tail1(1)));
        VssL = norm(vp1);
        
        vp2 = 1*Vp2(V2,dP2);
        v_tail2 = (v_tail0+vp2);
        V_tail2 = norm(v_tail2);
        Alpha_tail2 = rad2deg(atan2(v_tail2(3), v_tail2(1)));
        VssR = norm(vp2);
    end
    
    function [Z, X] = to_body_axes(L, D, Alpha)
        [R_V] = Wind_frame_matrix(Alpha);
        XYZ = R_V(1:3, 1:3)*[L;0;D];
        X = XYZ(1);
        Z = XYZ(3);
    end

end