function C = DC_cons_scen(x, dc, wind)
% C = AC_cons_scen(x, ac, wind, t_wind, j_des)
% returns an LMI with all the constraints
    
    C = [];
    
    % get indices for x
    PG_idx = 1:dc.N_G;
    Rus_idx = dc.N_G+1:2*dc.N_G;
    Rds_idx = 2*dc.N_G+1:3*dc.N_G;
    dus_idx = 3*dc.N_G+1:4*dc.N_G;
    dds_idx = 4*dc.N_G+1:5*dc.N_G;
    

    [N_t, N] = size(wind.P_m);
    
    % reshape x if it is a column vector
    if size(x, 2) == 1
        x = reshape(x, 5*dc.N_G, N_t);
    end
    
    for i = 1:N
        for t = 1:N_t
            
            % define reserve power
            R = x(dus_idx, t) * max(0, -wind.P_m(t, i)) ...
                            - x(dds_idx, t) * max(0, wind.P_m(t, i));                


            % define scenario power injection vector
            P_injs = dc.C_G * (x(PG_idx, t) + R) ...
                      + dc.C_w * wind.P_w(t, i) - dc.P_D(t, :)';

           
            % generator limits
            C = [C, dc.P_Gmin <= x(PG_idx, t) + R <= dc.P_Gmax];
            
            % bounds on R
            C = [C, -x(Rds_idx, t) <= R <= x(Rus_idx, t)];

            P_fs = dc.B_f * [dc.B_bustildeinv * P_injs(1:end-1); 0];

            % Line flow limits
            C = [C, -dc.P_fmax <= P_fs <= dc.P_fmax];
            
        end
    end
end