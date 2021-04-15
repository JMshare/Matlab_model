function [V, Alpha, Beta, Gamma] = Extract_aero_states(y)

    v = y(4:6,:);
    epsilon = y(10:12,:);
    
    V = sqrt(v(1,:).^2 + v(2,:).^2 + v(3,:).^2);
    alpha = atan2(v(3,:), v(1,:));
    Alpha = rad2deg(alpha);
    if V == 0
        beta = 0;
    else
        beta = -asin(v(2,:)./V);
    end
    Beta = rad2deg(beta);
    Gamma = rad2deg(epsilon(2,:)) - Alpha;
    
    
end