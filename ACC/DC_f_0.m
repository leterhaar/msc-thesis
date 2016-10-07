function C = DC_f_0(x, dc, wind, t_wind)
% returns the deterministic constraints

    PG_idx = 1:dc.N_G;
    Rus_idx = dc.N_G+1:2*dc.N_G;
    Rds_idx = 2*dc.N_G+1:3*dc.N_G;

    C = [];

    % define deterministic power injection vector  
    P_injf = (dc.C_G * x(PG_idx) - dc.P_D(t_wind, :)' + dc.C_w * wind.P_wf(t_wind));

    % power balance constraints
    C = [C, (ones(1, dc.N_b) * P_injf == 0):'Power bal f'];

    % generator limits
    C = [C, (dc.P_Gmin <= x(PG_idx) <= dc.P_Gmax):'Gen lims f'];

    % line flow limits
    C = [C,  (- dc.P_fmax <=  ...
               dc.B_f * [dc.B_bustildeinv * P_injf(1:end-1); 0] ...
               <= dc.P_fmax):'Line lims f'];

    % Non-negativity constraints for reserve requirements       
    C = [C, (x(Rus_idx) >= 0):'NN Rup', (x(Rds_idx) >= 0):'NN Rdown'];
end