%% Main function to generate tests
% type "runtests" to run
function tests = Test_functions
tests = functiontests(localfunctions);
end

%% Test Functions
function test_Rot_and_Euler(testCase)
    % Tests R and E matrices and their inverses based on exact expressions from Etkin book
    clear all;
    e = [0.1,0.2,0.3];
    [R_inv, E_inv, ~, ~, R, E, ~, ~] = Rotation_and_Euler_Matrices(e);
    
    R_ex = [cos(e(2))*cos(e(3)), cos(e(2))*sin(e(3)), -sin(e(2));
            -cos(e(1))*sin(e(3)) + sin(e(1))*sin(e(2))*cos(e(3)), cos(e(1))*cos(e(3)) + sin(e(1))*sin(e(2))*sin(e(3)), sin(e(1))*cos(e(2));
            sin(e(1))*sin(e(3)) + cos(e(1))*sin(e(2))*cos(e(3)), -sin(e(1))*cos(e(3)) + cos(e(1))*sin(e(2))*sin(e(3)), cos(e(1))*cos(e(2))];
    
    E_ex = [1, 0, -sin(e(2));
            0, cos(e(1)), sin(e(1))*cos(e(2));
            0, -sin(e(1)), cos(e(1))*cos(e(2))];
        
    E_ex_inv = [1, sin(e(1))*tan(e(2)), cos(e(1))*tan(e(2));
                0, cos(e(1)), -sin(e(1));
                0, sin(e(1))/cos(e(2)), cos(e(1))/cos(e(2))];
    
    assert( norm(R - R_ex)  < 1e-8 , 'R doesn''t agree with Etkin');
    assert( norm(R_inv - R_ex') < 1e-8 , 'R_inv doesn''t agree with Etkin');
    assert( norm(E - E_ex) < 1e-8 , 'E doesn''t agree with Etkin');
    assert( norm(E_inv - E_ex_inv)< 1e-8 , 'E_inv doesn''t agree with Etkin');
    
end

function test_dEuler(testCase)
    % Tests derivative of E and E_inv based on analytical derivation of the
    % expressions from Etkin book.
    clear all;
    e = [0.1,0.2,0.3];
    [~, ~, ~, dE_inv, ~, ~, ~, dE] = Rotation_and_Euler_Matrices(e);
    
    dE_ex = [[0, 0, 0; 0, -sin(e(1)), cos(e(1))*cos(e(2)); 0, -cos(e(1)), -sin(e(1))*cos(e(2))];
             [0, 0, -cos(e(2)); 0, 0, -sin(e(1))*sin(e(2)); 0, 0, -cos(e(1))*sin(e(2))];
             [0, 0, 0; 0, 0, 0; 0, 0, 0]];
         
    dE_inv_ex = [[0, cos(e(1))*tan(e(2)), -sin(e(1))*tan(e(2)); 0, -sin(e(1)), -cos(e(1)); 0, cos(e(1))/cos(e(2)), -sin(e(1))/cos(e(2))];
                 [0, sin(e(1))/(cos(e(2))^2), cos(e(1))/(cos(e(2))^2); 0, 0, 0; 0, sin(e(1))*(sin(e(2))/(cos(e(2))^2)), cos(e(1))*(sin(e(2))/(cos(e(2))^2))];
                 [0, 0, 0; 0, 0, 0; 0, 0, 0]];
    
    assert( norm(dE - dE_ex) < 1e-8 , 'dE doesn''t agree with Etkin');
    assert( norm(dE_inv - dE_inv_ex)< 1e-8 , 'dE_inv doesn''t agree with Etkin');

end

function test_dRot_and_dEuler_num(testCase)
    % Tests derivatives of the E and R and their inverses based on
    % numerical differentiation
    clear all;
    
    e = [0.1,0.2,0.3];
    dx = 1e-5;
    n = length(e);
    
    % dR_inv
    dR_inv = zeros(0,n);
    for i = 1:n
        e_pl = e;
        e_mn = e;
        e_pl(i) = e_pl(i)+dx;
        e_mn(i) = e_mn(i)-dx;
        [R_inv_pl, ~, ~, ~, ~, ~] = Rotation_and_Euler_Matrices(e_pl);
        [R_inv_mn, ~, ~, ~, ~, ~] = Rotation_and_Euler_Matrices(e_mn);
        diff = (R_inv_pl-R_inv_mn)/(2*dx);
        dR_inv = [dR_inv; diff];
    end
    
    [~, ~, dR_inv_an, ~, ~, ~] = Rotation_and_Euler_Matrices(e);
    assert( norm(dR_inv - dR_inv_an) < 1e-8 , 'dR_inv doesn''t agree numeric');
    
    
    % dE_inv
    dE_inv = zeros(0,n);
    for i = 1:n
        e_pl = e;
        e_mn = e;
        e_pl(i) = e_pl(i)+dx;
        e_mn(i) = e_mn(i)-dx;
        [~, E_inv_pl, ~, ~, ~, ~] = Rotation_and_Euler_Matrices(e_pl);
        [~, E_inv_mn, ~, ~, ~, ~] = Rotation_and_Euler_Matrices(e_mn);
        diff = (E_inv_pl-E_inv_mn)/(2*dx);
        dE_inv = [dE_inv; diff];
    end
    
    [~, ~, ~, dE_inv_an, ~, ~] = Rotation_and_Euler_Matrices(e);
    tot = (dE_inv_an+dE_inv)/2;
    tot(abs(tot)<1e-10) = 0;
    bin = (dE_inv_an - dE_inv)./tot;
    bin(isnan(bin)) = 0;
    bin(isinf(bin)) = 0;
    assert( norm(bin) < 1e-9 , 'dE_inv doesn''t agree numeric');
    %assert( norm(dE_inv - dE_inv_an) < 1e-8 , 'dE_inv doesn''t agree numeric');
    
    
    % dR
    dR = zeros(0,n);
    for i = 1:n
        e_pl = e;
        e_mn = e;
        e_pl(i) = e_pl(i)+dx;
        e_mn(i) = e_mn(i)-dx;
        [~, ~, ~, ~, R_pl, ~] = Rotation_and_Euler_Matrices(e_pl);
        [~, ~, ~, ~, R_mn, ~] = Rotation_and_Euler_Matrices(e_mn);
        diff = (R_pl-R_mn)/(2*dx);
        dR = [dR; diff];
    end
    
    [~, ~, ~, ~, ~, ~, dR_an] = Rotation_and_Euler_Matrices(e);
    assert( norm(dR - dR_an) < 1e-8 , 'dR doesn''t agree numeric');
    
    
    % dE
    dE = zeros(0,n);
    for i = 1:n
        e_pl = e;
        e_mn = e;
        e_pl(i) = e_pl(i)+dx;
        e_mn(i) = e_mn(i)-dx;
        [~, ~, ~, ~, ~, E_pl] = Rotation_and_Euler_Matrices(e_pl);
        [~, ~, ~, ~, ~, E_mn] = Rotation_and_Euler_Matrices(e_mn);
        diff = (E_pl-E_mn)/(2*dx);
        dE = [dE; diff];
    end
    
    [~, ~, ~, ~, ~, ~, ~, dE_an] = Rotation_and_Euler_Matrices(e);
    assert( norm(dE - dE_an) < 1e-8 , 'dE doesn''t agree numeric');
    
end

% the following usage of Forces_and_Moments and global y0 u0 allocaton is
% depreciated
% function test_dF_num(testCase)
%     % tests dF numerically
%     clear all;
%     Load_params(1);
%     global y0 u0
%     y = y0;
%     u = u0;
%     [~, ~, dF_an] = Forces_and_Moments(y,u); 
%     
%     dx = 1e-5;
%     n = length(y);
%     dF = 0*dF_an;
%     for i = 1:n
%         y_pl = y;
%         y_mn = y;
%         y_pl(i) = y_pl(i)+dx;
%         y_mn(i) = y_mn(i)-dx;
%         F_pl = Forces_and_Moments(y_pl,u);
%         F_mn = Forces_and_Moments(y_mn,u);
%         diff = (F_pl-F_mn)/(2*dx);
%         dF(:,i) = diff;
%     end
%     
%     tot = (dF_an+dF)/2;
%     tot(abs(tot)<1e-10) = 0;
%     bin = (dF_an - dF)./tot;
%     bin(isnan(bin)) = 0;
%     bin(isinf(bin)) = 0;
%     assert( norm(bin) < 1e-7 , 'dF doesn''t agree numeric');
%     
% end


function test_cross_matrix_num(testCase)
    % tests the cross product matrix derivatives numerically
    clear all;
    omega = [0.03, 0.05, 0.08]';
    [~, dOmg] = Cross_Matrix(omega);
    
    dx = 1e-5;
    n = length(omega);
    dOmg_num = zeros(0,n);
    for i = 1:n
        omega_pl = omega;
        omega_mn = omega;
        omega_pl(i) = omega_pl(i)+dx;
        omega_mn(i) = omega_mn(i)-dx;
        diff = (Cross_Matrix(omega_pl)-Cross_Matrix(omega_mn))/(2*dx);
        dOmg_num = [dOmg_num; diff];
    end
    
    bin = (dOmg - dOmg_num)./((dOmg+dOmg_num)/2);
    bin(isnan(bin)) = 0;
    assert( norm(bin) < 1e-9 , 'dOmg doesn''t agree numeric');
end



