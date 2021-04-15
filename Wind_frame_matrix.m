function [R_V] = Wind_frame_matrix(Alpha)
    % Transformation matrix to and from body frame
    % expects angle of attack in degs
    R_V = eye(6);
    R_V(1,1) = +sind(Alpha);
    R_V(1,3) = -cosd(Alpha);
    R_V(3,1) = -cosd(Alpha);
    R_V(3,3) = -sind(Alpha);
end