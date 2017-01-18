function cost = a_posteriori_cost(P_G, R)
% cost = a_posteriori_cost(P_G, R)
% 
% Calculates the a posteriori cost, based on the values for P_G and R
% Gets ac, wind, t from the caller

    check_class({P_G, R}, {'double', 'double'});
    
    ac = evalin('caller', 'ac');
    
    % deterministic cost
    cost = ac.c_li'*P_G + P_G' * diag(ac.c_qu) * P_G;
    
    % determine maximum up and downspinning reserve
    total_reserves = sum(R, 1);
    R_us = R(:, arg_max(total_reserves));
    R_ds = R(:, arg_min(total_reserves));
    
    % add reserve cost
    cost = cost + ac.c_ds' * R_ds + ac.c_us' * R_us;
end
