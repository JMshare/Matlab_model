function [dy] = ODE_fun(y,u)
% ODE function dy = f(y,u) representing the nonlinear equations of motion for
% a 3D airplane.
    
    global m I minv Iinv c0
    
    % Extract the state variable groups
    x = y(1:3); % [x, y, z] in the inertial frame
    v = y(4:6); % [u, v, w] in the body frame
    omega = y(7:9); % [p, q, r] in the body frame
    epsilon = y(10:12); % [phi, theta, psi] in the 'euler' frame

    % Create a cross product matrix
    [Omg] = Cross_Matrix(omega);
    
    % Create the transformation matrices
    [Rinv, Einv] = Rotation_and_Euler_Matrices(epsilon);
      
    % Create the forcing terms vectors
    [F, M, hp] = Forces_and_Moments(y,u);
    if (norm(hp) > 0) && (norm(omega) > 0)
        []; % debug point
    end

    % Equations of motions for the individual state groups
    dx = Rinv*(v + c0);
    dv = minv*(F - Omg*m*v);
    domega = Iinv*(M - Omg*(I*omega+hp));
    depsilon = Einv*omega;
    

    % Equations of motions vector
    dy = [dx;
          dv;
          domega;
          depsilon];

end


