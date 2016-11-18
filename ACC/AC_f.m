function Obj = AC_f(x, ac, wind, t)
% returns the objective function

    % extract  variables from x
    W_f = x{1};
    R = x{4};

    % define indexes for variables
    Rus_idx = 1:ac.N_G;
    Rds_idx = ac.N_G+1:2*ac.N_G;
    
    Obj = 0;
    lambda = 1;

    % add power generation costs
    for j = 1:ac.N_G
        k = ac.Gens(j);
        Obj = Obj + ac.c_us(j)*(trace(ac.Y_k(k)*W_f) + ac.P_D(t,k) ...
                                                -ac.C_w(k)*wind.P_wf(t));
    end

    % add reserve requirements costs
    Obj = Obj + lambda * (ac.c_us' * R(Rus_idx)  ...
                       + ac.c_ds' * R(Rds_idx));
end