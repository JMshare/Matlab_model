function [yn, un, cstar] = Trim_conditions(Vstar, cslast)
%%% Trim the airplane %%%
    
    % Specify conditions
    cs0 = [Vstar/2; Vstar/2; 0; 0.5; 0.5]; % [u, w; dE, dP1=dP2, dP3]
    cs1 = [Vstar/2; Vstar/2; 0; 1; 0];
    cs2 = [Vstar/2; Vstar/2; 0; 0; 1];
    
    cs = [cs0, cs1, cs2, cslast];
  
    lb = [0; -Vstar; -1; 0; 0];
    ub = [Vstar; Vstar; 1; 1; 1];

    uAS = 0;
    
    % Optimize Trim
    options = optimoptions('fmincon', 'StepTolerance', 1e-16, 'MaxFunctionEvaluations', 1000);
    Js = 1e12*ones(size(cs,2),1);
    cstar = 0*cs;
    for i = 1:length(Js)
        c = cs(:,i);
        cstar(:,i) = fmincon(@Trim_obj_fun, c, [], [], [], [], lb, ub, @nlcon, options);
        Js(i) = Trim_obj_fun(cstar(:,i));
        if Js(i) < 1e-1
            break;
        end
    end
    [~, idx] = min(Js);
    cstar = cstar(:,idx);
    [yn, un] = insert_c(cstar);
    [Fn, Mn] = Forces_and_Moments(yn,un);
    [Fn([1,3]);Mn(2)]

    %% Function definitions
    function [J] = Trim_obj_fun(c)
        % Requires zero forces and moments on trim
        [y, u] = insert_c(c);
        [F, M] = Forces_and_Moments(y,u);
        FM = [F([1,3]); M(2)]; % only the longit trim
        Lam = 1e-2*diag([0,0,0,0,0]);
        % J = FM'*FM + c'*Lam*c + 1*c(4)*c(5) + 0*c(4) +
        % 0*1e-2*c(4)*(u(1)>20) + atan2(c(2), c(1))*(u(1)<20); really nice
        % result but hard to justify this
        J = FM'*FM + c'*Lam*c + 1*c(4)*c(5) + 2*c(4) + 2*c(5) + 0*1e-2*c(4)*(u(1)>20) + atan2(c(2), c(1)) + 10*c(4)*(c(1)>14.5); % this does the job even though I may wana tune the weights a bit to get the transitions where i need them
    end

    function [y, u] = insert_c(c)
        y = [0; 0; 0; c(1); 0; c(2); 0; 0; 0; 0; atan2(c(2), c(1)); 0]; % u and w
        u = [0; c(3); 0; uAS; 0; c(4); c(4); c(5); 0]; % props
    end

    function [hc, gc] = nlcon(c)
        hc = [rad2deg(atan2(c(2), c(1))) - 91; 
             -rad2deg(atan2(c(2), c(1))) - 8;
             (sqrt(c(1)^2 + c(2)^2)*sign(c(1)) - Vstar) - 0.2;
             -(sqrt(c(1)^2 + c(2)^2)*sign(c(1)) - Vstar + 0.2)]; % inequality constraint h(c) <= 0
        gc = 0*[sqrt(c(1)^2 + c(2)^2)*sign(c(1)) - Vstar]; % equality constraint g(c) == 0
    end

end
