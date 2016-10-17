function Obj = DC_f_obj(x, dc, wind, t_wind)
% returns the objective function

    % define indexes for variables
    PG_idx = 1:dc.N_G;
    Rus_idx = dc.N_G+1:2*dc.N_G;
    Rds_idx = 2*dc.N_G+1:3*dc.N_G;
    
    Obj = 0;
    lambda = 1;

    % add power generation costs
    Obj = Obj + dc.c_li' * x(PG_idx) ...
                        + x(PG_idx)' * diag(dc.c_qu) * x(PG_idx);

    % add reserve requirements costs
    Obj = Obj + lambda * (dc.c_us' * x(Rus_idx)  ...
                       + dc.c_ds' * x(Rds_idx));
end