function C = feasibleW(W_f, W_s, R_us, R_ds, Pw)
% C = feasible_set(W, Pw)
% returns a set of constraints that determine the feasible set for a state
% of the network
% 
% ARGUMENTS
% =========
% W : sdpvar with some state of the network
% R : sdpvar with [Rds; Rus]
% Pw : vector with some wind realization

    % get network and wind model and wind time
    ac = evalin('caller', 'ac');
    wind = evalin('caler', 'wind');
    t = evalin('caller', 't');
    
    C = [];
    % define constraints for every bus
    for k = 1:ac.N_b
        % P_inj (1)
        C = [C, ac.P_min(k) ...
            <= trace(ac.Y_k(k)*W_s) + ac.P_D(t, k) - ac.C_w(k)*Pw(t) <= ...
               ac.P_max(k)];

        % Q_inj (2)
        C = [C, ac.Q_min(k) ...
            <= trace(ac.Ybar_k(k)*W_s) + ac.Q_D(t, k) <=...
               ac.Q_max(k)];

        % V_bus (3)
        C = [C, ac.V_min(k)^2 ...
            <= trace(ac.M_k(k)*W_s) <= ...
               ac.V_max(k)^2];
    end
    
    
    for j = 1:ac.N_G
        k = ac.Gens(j);
        % Bound R between R_us and R_ds
        C = [C, -R_ds(j) ...
             <= trace(ac.Y_k(k)*(W_s-W_f)) - ac.C_w(k)*(P_w(t) - wind.P_wf(t)) <= ...
                R_us(j+ac.N_G)];

    end
end