function C = DC_cons_det(x, dc, wind)
% DC_cons_det(x, dc, wind)
% returns the deterministic constraints

    PG_idx = 1:dc.N_G;
    Rus_idx = dc.N_G+1:2*dc.N_G;
    Rds_idx = 2*dc.N_G+1:3*dc.N_G;
    dus_idx = 3*dc.N_G+1:4*dc.N_G;
    dds_idx = 4*dc.N_G+1:5*dc.N_G;

    C = [];
    N_t = size(wind.P_m, 1);
    
    % reshape x if it is a column vector
    if size(x, 2) == 1
        x = reshape(x, 5*dc.N_G, N_t);
    end
    
    for t = 1:N_t

        % define deterministic power injection vector  
        P_injf = (dc.C_G * x(PG_idx, t) - dc.P_D(t, :)' ...
                 + dc.C_w * wind.P_wf(t));

        % power balance constraints
        C = [C, (ones(1, dc.N_b) * P_injf == 0):...
                                        sprintf('Power bal det t%2i', t)];

                                    
        % generator limits
        C = [C, (dc.P_Gmin <= x(PG_idx, t) <= dc.P_Gmax):...
                                         sprintf('Gen lims det t%2i', t)];

        % line flow limits
        C = [C,  (- dc.P_fmax <=  ...
                   dc.B_f * [dc.B_bustildeinv * P_injf(1:end-1); 0] ...
                   <= dc.P_fmax):sprintf('Line lims det t%2i', t)];

        % Non-negativity constraints for reserve requirements       
        C = [C, (x(Rus_idx, t) >= 0):sprintf('NN Rup t%2i', t), ...
                (x(Rds_idx, t) >= 0):sprintf('NN Rdown t%2i', t)];

        C = [C, (ones(1,dc.N_G)*x(dus_idx, t) == 1):sprintf('Balu t%2i',t), ...
                (ones(1,dc.N_G)*x(dds_idx, t) == 1):sprintf('Bald t%2i',t)];
            
        
    end
end