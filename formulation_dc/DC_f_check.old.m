function [params_act, residuals] = DC_f_check(x, i, dc, wind, t_wind, j_des)
% [C; params] = DC_f_ineq(x, i, dc, wind, t_wind)
% Function to check a constraint for some scenario for solution x
    
    if nargin < 6
        j_des = 0;
    end

    tol = 1e-6;

    PG_idx = 1:dc.N_G;
    Rus_idx = dc.N_G+1:2*dc.N_G;
    Rds_idx = 2*dc.N_G+1:3*dc.N_G;
    dus_idx = 3*dc.N_G+1:4*dc.N_G;
    dds_idx = 4*dc.N_G+1:5*dc.N_G;
    
    % define reserve power
    R = x(dus_idx) * max(0, -wind.P_m(t_wind, i)) ...
                    - x(dds_idx) * max(0, wind.P_m(t_wind, i));                
                
                
    % define scenario power injection vector
    P_injs = dc.C_G * (x(PG_idx) + R) ...
              + dc.C_w * wind.P_w(t_wind, i) - dc.P_D(t_wind, :)';
    N_j = 4*dc.N_G + 2*dc.N_l;
    
    % preallocate residuals and parameters
    params = zeros(N_j, 2);
    residuals = zeros(N_j, 1);
    
    j = 1;

    for k = 1:dc.N_G   

        % generator lower limit
        if j_des == j || j_des == 0
            residuals(j) = - dc.P_Gmin(k) + x(PG_idx(k)) + R(k);
            params(j, :) = [i j];
        end
        j = j+1;

        % generator upper limit
        if j_des == j || j_des == 0
            residuals(j) = - x(PG_idx(k)) - R(k) + dc.P_Gmax(k);
            params(j, :) = [i j];
        end
        j = j+1;    

        % R lower bound
        if j_des == j || j_des == 0
            residuals(j) = x(Rds_idx(k)) + R(k);
            params(j, :) = [i j];
        end
        j = j+1;   

        % R upper bound
        if j_des == j || j_des == 0
            residuals(j) = - R(k) + x(Rus_idx(k));
            params(j, :) = [i j];
        end
        j = j+1;   

    end
    
    P_fs = dc.B_f * [dc.B_bustildeinv * P_injs(1:end-1); 0];
    
    for k = 1:dc.N_l
        
        % line flow lower limit
        if j_des == j || j_des == 0
            residuals(j) = dc.P_fmax(k) + P_fs(k);
            params(j, :) = [i j];
        end
        j = j+1;   
        
        % line flow lower limit
        if j_des == j || j_des == 0
            residuals(j) = - P_fs(k) + dc.P_fmax(k);
            params(j, :) = [i j];
        end
        j = j+1;   
        
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

    
    

