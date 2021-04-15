function [A_num, B_num, fAero] = ODE_Jac_num(y, u)
    %% Numerically determine the Jacobian of the ODE function ydot = f(y,u)
    % at state y and control u.
    
    dx = 1e-5;
    n = length(y);
    m = length(ODE_fun(y,u));
    A_num = zeros(m,n);
    for i = 1:n
        if i==12
            [];
        end
        y_pl = y;
        y_mn = y;
        y_pl(i) = y_pl(i) + dx;
        y_mn(i) = y_mn(i) - dx;
        diff = (ODE_fun(y_pl,u) - ODE_fun(y_mn,u))/(2*dx);
        A_num(:,i) = diff;
    end
    
    dx = 1e-3;
    n = length(u);
    m = length(ODE_fun(y,u));
    B_num = zeros(m,n);
    for i = 1:n
        if i == 5
            [];
        end
        if i >= 6
            dx = 0.1;
        else
            dx = 1e-3;
        end
        u_pl = u;
        u_mn = u;
        u_pl(i) = u_pl(i) + dx;
        u_mn(i) = u_mn(i) - dx;
        diff = (ODE_fun(y,u_pl) - ODE_fun(y,u_mn))/(2*dx);
        B_num(:,i) = diff;
    end

    
    c = u;
    dx = 1e-5;
    [V, Alpha, Beta] = Extract_aero_states(y);
    Vpl = V + dx;
    Vmn = V - dx;
    [u,v,w] = Extract_body_states(Vpl, Alpha, Beta);
    y_pl = inputate(u,v,w,y);
    [u,v,w] = Extract_body_states(Vmn, Alpha, Beta);
    y_mn = inputate(u,v,w,y);
    fV = (ODE_fun(y_pl,c) - ODE_fun(y_mn,c))/(2*dx);

    dx = 1e-5;
    [V, Alpha, Beta] = Extract_aero_states(y);
    Alphapl = deg2rad(Alpha) + dx;
    Alphamn = deg2rad(Alpha) - dx;
    [u,v,w] = Extract_body_states(V, rad2deg(Alphapl), Beta);
    y_pl = inputate(u,v,w,y);
    [u,v,w] = Extract_body_states(V, rad2deg(Alphamn), Beta);
    y_mn = inputate(u,v,w,y);
    fAlpha = (ODE_fun(y_pl,c) - ODE_fun(y_mn,c))/(2*dx);

    dx = 1e-5;
    [V, Alpha, Beta] = Extract_aero_states(y);
    Betapl = deg2rad(Beta) + dx;
    Betamn = deg2rad(Beta) - dx;
    [~,v,~] = Extract_body_states(V, Alpha, rad2deg(Betapl));
    y_pl = inputate(y(4),v,y(6),y);
    [~,v,~] = Extract_body_states(V, Alpha, rad2deg(Betamn));
    y_mn = inputate(y(4),v,y(6),y);
    fBeta = (ODE_fun(y_pl,c) - ODE_fun(y_mn,c))/(2*dx);

    fAero = [fV, fAlpha, fBeta];
    
end

function [y] = inputate(u,v,w,y)
    y(4) = u;
    y(5) = v;
    y(6) = w;
end
