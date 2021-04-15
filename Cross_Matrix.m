function [Omg, dOmg] = Cross_Matrix(omega)
% Returns the cross product matrix Omg to represent (omega x v) as (Omg*v)
% Also returns its derivative wrt omega, being 3 matrices stored as
% [A;B;C].

    Omg = [ 0        -omega(3)  omega(2);
            omega(3)  0        -omega(1);
           -omega(2)  omega(1)  0       ];
   
    
    % Derivatives
    dOmgx = [0, 0, 0;
             0, 0, -1;
             0, 1, 0];
    dOmgy = [0, 0, 1;
             0, 0, 0;
             -1, 0, 0];
    dOmgz = [0, -1, 0;
             1, 0, 0;
             0, 0, 0];
    dOmg = [dOmgx; 
            dOmgy; 
            dOmgz];
end