function C = feasibleW(W_s, P_w)
% C = feasible_set(W, Pw)
% returns a set of constraints that determine the feasible set for a state
% of the network
% 
% ARGUMENTS
% =========
% W     : sdpvar with some state of the network
% Pw    : vector with some wind realization
%   
% RETURNS
% =======
% C     : the feasible set

    check_class({W_s, P_w}, {'sdpvar', 'double|sdpvar'});

    % get network and wind model and wind time
    ac = evalin('caller', 'ac');
    wind = evalin('caller', 'wind');
    t = evalin('caller', 't');
    
    % extend it to 24hours the same if only a scalar
    if length(P_w) == 1
        P_w = ones(24,1) * P_w;
    end
    
    C = [];
    % define constraints for every bus
    for k = 1:ac.N_b
        % P_inj (1)
        C = [C, ac.P_min(k) ...
            <= trace(ac.Y_k(k)*W_s) + ac.P_D(t, k) - ac.C_w(k)*P_w(t) <= ...
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
    
%     for l = 1:ac.N_l
%         % line constraints
%         C = [C, trace(ac.Y_lm(l) * W_s) <= ac.P_lmmax(l)];
%     end
    
    % refbus constraints
    refbus = ac.refbus + ac.N_b;    
    C = [C, W_s(refbus, refbus) == 0];
         
%     for j = 1:ac.N_G
%         k = ac.Gens(j);
%         % Bound R between R_us and R_ds
%         C = [C, -R_ds(j) ...
%              <= trace(ac.Y_k(k)*(W_s-W_f)) - ac.C_w(k)*(P_w(t) - wind.P_wf(t)) <= ...
%                 R_us(j)];
% 
%     end
end