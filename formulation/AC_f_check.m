function [params_act, residuals] = AC_f_check(x, i, ac, wind, t, j_des)
% [params, residuals] = AC_f_ineq(x, i, ac, wind, t_wind, [j_des])
% Function to check a constraint j_des for some scenario i for solution x
% j_des can also be left out, then all constraints for scenario i will
% be checked
% returns the parameters of the active constraints and all the residuals
    
    if nargin < 6
        j_des = 0;
    end

    tol = 1e-6;

    % extract variables from x
    W_f = x{1};
    W_m = x{2};
    R = x{3};
    R_us = R(1:ac.N_G);
    R_ds = R(ac.N_G+1:2*ac.N_G);
    
    dat dat zus zo
    N_j = 6*ac.N_b + 2*ac.N_G + 1;
    
    % preallocate residuals and parameters
    params = zeros(N_j, 2);
    residuals = zeros(N_j, 1);
    
    j = 1;
    
    for k = 1:ac.N_b
        
        % P_inj lower (1')
        if j == j_des || j_des == 0
            residuals(j) = -ac.P_min(k) ...
                + trace(ac.Y_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                         + ac.P_D(t, k) - ac.C_w(k)*wind.P_w(t, i);
            params(j,:) = [i j];
        end
        j = j + 1;
        
        % P_inj upper (1')
        if j == j_des || j_des == 0
            residuals(j) = -(trace(ac.Y_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                         + ac.P_D(t, k) - ac.C_w(k)*wind.P_w(t, i)) + ...
                   ac.P_max(k);
            params(j,:) = [i j];
        end
        j = j + 1;

        % Q_inj lower (2')
        if j == j_des || j_des == 0
            residuals(j) = -ac.Q_min(k) ...
                 + trace(ac.Ybar_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                            + ac.Q_D(t, k);
            params(j,:) = [i j];
        end
        j = j + 1;
        
        % Q_inj upper (2')
        if j == j_des || j_des == 0
            residuals(j) = -(trace(ac.Ybar_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                            + ac.Q_D(t, k)) + ac.Q_max(k);
            params(j,:) = [i j];
        end
        j = j + 1;


        % V_bus lower (3')
        if j == j_des || j_des == 0
            residuals(j) = -ac.V_min(k)^2 ...
                 + trace(ac.M_k(k)*(W_f + W_m * wind.P_m(t, i)));
            params(j,:) = [i j];
        end
        j = j + 1;   
        
        % V_bus upper (3')
        if j == j_des || j_des == 0
            residuals(j) = -trace(ac.M_k(k)*(W_f+W_m*wind.P_m(t, i))) + ...
                    ac.V_max(k)^2;
            params(j,:) = [i j];
        end
        j = j + 1;   

    end
        

    for g = 1:ac.N_G

        % bus index
        k = ac.Gens(g);

        % Lower bound R between R_us and R_ds
        if j == j_des || j_des == 0
            residuals(j) = R_ds(g) ...
                 + trace(ac.Y_k(k)*(W_m*wind.P_m(t, i))) ...
                                    - ac.C_w(k)*wind.P_m(t, i);
            params(j,:) = [i j];
        end
        j = j + 1;
        
        % Lower bound R between R_us and R_ds
        if j == j_des || j_des == 0
            residuals(j) = -(trace(ac.Y_k(k)*(W_m*wind.P_m(t, i))) ...
                                    - ac.C_w(k)*wind.P_m(t, i)) + R_us(g);
            params(j,:) = [i j];
        end
        j = j + 1;
        
    end
    
    if j == j_des || j_des == 0
        residuals(j) = min(eig(W_f + W_m * wind.P_m(t, i)));
        params(j,:) = [i j];
    end
    
    if j_des == 0
        params_act = params((residuals < tol) & (residuals > -tol), :);
    else
        if residuals(j_des) < tol && residuals(j_des) > -tol
            params_act = params(j_des, :);
        else
            params_act = [];
        end
        residuals = residuals(j_des);
    end
end

