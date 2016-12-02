function [C, params] = DC_f_ineq(x, i, dc, wind, t_wind, j_des)
% [C, params] = DC_f_ineq(x, i, dc, wind, t_wind)
% Function to create a constraint set based on a scenario

    if nargin < 6
        j_des = 0;
    end

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
    C = [];
    params = zeros(4*dc.N_G + 2*dc.N_l, 2);
    j = 1;
   
    for k = 1:dc.N_G   

        % generator lower limit
        
        % TODO make all single sided
        if j == j_des || j_des == 0
            C = [C, (dc.P_Gmin(k) <= x(PG_idx(k)) + R(k))...
                    :sprintf('Gen lower s%i g%i', i, k)];
        end
        params(j, :) = [i j];
        j = j+1;

        % generator upper limit
        if j == j_des || j_des == 0
        C = [C, (x(PG_idx(k)) + R(k) <= dc.P_Gmax(k))...
                :sprintf('Gen upper s%i g%i', i, k)];
        end
        params(j, :) = [i j];
        j = j+1;    

        % R lower bound
        if j == j_des || j_des == 0
        C = [C, (-x(Rds_idx(k)) <= R(k)):sprintf('R lower s%i g%i',i, k)];
        end
        params(j, :) = [i j];
        j = j+1;   

        % R upper bound
        if j == j_des || j_des == 0
        C = [C, (R(k) <= x(Rus_idx(k))):sprintf('R upper s%i g%i',i, k)];
        end
        params(j, :) = [i j];
        j = j+1;   

    end
    
    
    P_fs = dc.B_f * [dc.B_bustildeinv * P_injs(1:end-1); 0];
    
    for k = 1:dc.N_l
        
        % line flow lower limit
        if j == j_des || j_des == 0
        C = [C, (- dc.P_fmax(k) <= P_fs(k)):sprintf('Line lower s%i l%i',i,k)];
        end
        params(j, :) = [i j];
        j = j+1;   
        
        % line flow lower limit
        if j == j_des || j_des == 0
        C = [C, (P_fs(k) <= dc.P_fmax(k)):sprintf('Line upper s%i l%i',i,k)];
        end
        params(j, :) = [i j];
        j = j+1;   

    end
    
end
