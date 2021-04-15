function [R_inv, E_inv, dR_inv, dE_inv, R, E, dR, dE] = Rotation_and_Euler_Matrices(epsilon)
% Returns the rotational matrix of Euler rotations from original to rotated
% position, 
% and a transformation matrix E to transform de in the 'euler' frame 
% to omega in the body frame.
%
% R is a matrix of rotations therefore inv(R) = R'.
% This doesn't hold for E, in fact E can be even singular (gimbal lock).
%
%
% Oppositely, looking at a vector x in the inertial frame with
% origin at the center or rotation pointing to a point mass in the body,
% the transformation of this vector due to the body rotation can be
% represented as x(1) = R*x(0). 
% We can also find the omega from the cross product matrix Omg by 
% expressing v(1) by differentiation as 
% v(1) = dR*x(0) = dR*R'*x(0) = Omg*x(t) = (omega cross x).
%
% Alternatvely, we can express the omega in the body frame directly from
% the Eulers rotations as omega = E*depsilon
%
% Input the euler angles e = epsilon = [phi, theta, psi] in [rad]
%

    % Prepare some blocks
    Rx = @(a)   [1         0         0;
                 0       cos(a) -sin(a);
                 0       sin(a)  cos(a)];
    Ry = @(a)   [cos(a)  0 sin(a);
                 0       1 0     ;
                 -sin(a) 0 cos(a)];
    Rz = @(a)   [cos(a) -sin(a) 0;
                 sin(a)  cos(a) 0;
                 0         0    1];
             
    Rxp = @(a) [0     0          0    ;
                0     -sin(a)  -cos(a);
                0     cos(a)   -sin(a)];
    Ryp = @(a) [-sin(a)  0  cos(a);
                0        0  0       ;
                -cos(a)  0  -sin(a)];
    Rzp = @(a) [-sin(a) -cos(a) 0;
                cos(a)  -sin(a) 0;
                0          0    0];
            
    Rxpp = @(a) [0     0          0    ;
                0     -cos(a)  sin(a);
                0     -sin(a)   -cos(a)];
    Rypp = @(a) [-cos(a)  0  -sin(a);
                0        0  0       ;
                sin(a)  0  -cos(a)];
    Rzpp = @(a) [-cos(a) sin(a) 0;
                -sin(a)  -cos(a) 0;
                0          0    0];
    
    
        
    x = [1; 0; 0];
    y = [0; 1; 0];
    z = [0; 0; 1];
    e = epsilon; % (phi, theta, psi)

    % Assume the Z-Y-X rotation (Yaw, Pitch, Roll) by (psi, theta, phi)
    R_inv = Rz(e(3))*Ry(e(2))*Rx(e(1));
    
    R = Rx(e(1))'*Ry(e(2))'*Rz(e(3))';
    
    E = [x, Rx(e(1))'*y, Rx(e(1))'*Ry(e(2))'*z]; 
    E_inv = pinv(E);
    
    
    % Derivatives
    dR_inv = [Rz(e(3))*Ry(e(2))*Rxp(e(1)); 
              Rz(e(3))*Ryp(e(2))*Rx(e(1)); 
              Rzp(e(3))*Ry(e(2))*Rx(e(1))];
          
    
    dR = [Rxp(e(1))'*Ry(e(2))'*Rz(e(3))';
          Rx(e(1))'*Ryp(e(2))'*Rz(e(3))';
          Rx(e(1))'*Ry(e(2))'*Rzp(e(3))'];
          
    
    dE = [[0*x, Rxp(e(1))'*y, Rxp(e(1))'*Ry(e(2))'*z];
          [0*x, 0*y,          Rx(e(1))'*Ryp(e(2))'*z];
          [0*x, 0*y,          0*z]];
      
    dinv = pinv(E);
    dE_inv = [-dinv*[0*x, Rxp(e(1))'*y, Rxp(e(1))'*Ry(e(2))'*z]*dinv;
              -dinv*[0*x, 0*y,          Rx(e(1))'*Ryp(e(2))'*z]*dinv;
              -dinv*[0*x, 0*y,          0*z]*dinv];
      
end