function Obj = AC_f(x, ac, wind, t)
% returns the objective function

    % extract  variables from x
    W_f = x{1};
    W_mus = x{2};
    W_mds = x{3};
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
%         for i = [argmax(wind.P_m(t, :)) arg_min(wind.P_m(t, :))]
%             
%             % define scenario state
%             W_s = W_f + W_mus * max(0, -wind.P_m(t, i)) ...
%                       - W_mds * max(0, wind.P_m(t, i));
%                   
%             Obj = Obj + ac.c_us(j)*trace(ac.Y_k(k)*W_s-W_f) - ac.C_w(k)*wind.P_m(t, i);
%         end
    end
    
    

    % add reserve requirements costs
    Obj = Obj + lambda * (ac.c_us' * R(Rus_idx)  ...
                       + ac.c_ds' * R(Rds_idx));

    
end