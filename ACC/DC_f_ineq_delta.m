function C = DC_f_ineq_delta(x, delta, dc, t_wind)
% C = DC_f_ineq_delta(x, delta, dc, t_wind)
%
% delta(1) = Pw, delta(2) = max(0, -Pm), delta(3) = max(0, Pm)
% Function to create a constraint set based on a scenario

    PG_idx = 1:dc.N_G;
    Rus_idx = dc.N_G+1:2*dc.N_G;
    Rds_idx = 2*dc.N_G+1:3*dc.N_G;
    dus_idx = 3*dc.N_G+1:4*dc.N_G;
    dds_idx = 4*dc.N_G+1:5*dc.N_G;
    
    % define reserve power
    R = x(dus_idx) * delta(2) - x(dds_idx) * delta(3);                
                
                
    % define scenario power injection vector
    P_injs = dc.C_G * (x(PG_idx) + R) ...
              + dc.C_w * delta(1) - dc.P_D(t_wind, :)';
    C = [];
    params = zeros(4*dc.N_G + 2*dc.N_l, 2);
    j = 1;
    i = 1;
    j_des = 0;
    for k = 1:dc.N_G   

        % generator lower limit
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
