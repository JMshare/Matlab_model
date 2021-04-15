function [u,v,w] = Extract_body_states(V, Alpha, Beta)

    u = V*cosd(Alpha); % only when v==0
    w = V*sind(Alpha); % only when v==0
    v = -sind(Beta)*V;

end