function [uN] = Trim_conditions_later(y,u)
%%% Lateral trim the airplane %%%
    
    % Specify conditions
    c0 = [0;0]; % [A, R]
  
    lb = [-1; -1];
    ub = [ 1;  1];
    
    % Optimize Trim
    options = optimoptions('fmincon', 'StepTolerance', 1e-16, 'MaxFunctionEvaluations', 1000, 'OptimalityTolerance', 1e-16);
    cstar = fmincon(@Trim_obj_fun, c0, [], [], [], [], lb, ub, [], options);
    [uN] = insert_c(cstar);

    %% Function definitions
    function [J] = Trim_obj_fun(c)
        % Requires zero forces and moments on trim
        [u2] = insert_c(c);
        [F, M] = Forces_and_Moments(y,u2);
        FM = [F([]); M([1,3])]; % only the later trim
        J = FM'*FM + 0.01*(c'*c) + 0*c(2)*c(2); 
    end

    function [u2] = insert_c(c)
        u2 = u + [c(1); 0; c(2); 0; 0.0*c(1); -0.0*c(2); 0.0*c(2); 0; 0]; 
    end

end
