function C = DC_f_eq(x, i, dc, wind, t_wind)
% C = DC_f_eq(x, dc, wind, t_wind)
% returns the equality constraints for a certain scenario

    PG_idx = 1:dc.N_G;
    Rus_idx = dc.N_G+1:2*dc.N_G;
    Rds_idx = 2*dc.N_G+1:3*dc.N_G;
    dus_idx = 3*dc.N_G+1:4*dc.N_G;
    dds_idx = 4*dc.N_G+1:5*dc.N_G;
    
    % power balance constraints in for scenario i
    R = x(dus_idx) * max(0, -wind.P_m(t_wind, i)) ...
                    - x(dds_idx) * max(0, wind.P_m(t_wind, i));
    P_injs = dc.C_G * (x(PG_idx) + R) ...
          + dc.C_w * wind.P_w(t_wind, i) - dc.P_D(t_wind, :)';

    C = (ones(1, dc.N_b) * P_injs == 0):sprintf('Power bal s%i', i);
    
    
end