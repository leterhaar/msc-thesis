function C = feasibleW(W_s, P_w, ac, notation)
% C = feasible_set(W, Pw)
% returns a set of constraints that determine the feasible set for a state
% of the network
% 
% ARGUMENTS
% =========
% W     : sdpvar with some state of the network
% Pw    : vector with some wind realization
% ac    : instance of AC_model or subnetwork 
%         (optional, if left empty 'ac' will be evaluated from caller)
% notation : either 'rectangular' or 'complex', optional, default
%            'rectangular'
%   
% RETURNS
% =======
% C     : the feasible set

    if nargin < 4
        notation = 'rectangular';
    end

    check_class({W_s, P_w}, {'sdpvar', 'double|sdpvar'});

    % get network and wind model and wind time
    if nargin < 3
        ac = evalin('caller', 'ac');
    end
    t = evalin('caller', 't');
    
    % extend it to 24hours the same if only a scalar
    if length(P_w) == 1
        P_w = ones(24,1) * P_w;
    end
    
    C = [];
    
    if strcmpi(notation, 'rectangular')
        % define constraints for every bus
        for k = 1:ac.N_b
            
            if not(isnan(ac.P_min(k)))
                
                % P_inj (1)
                C = [C, (ac.P_min(k) ...
                    <= trace(ac.Y_k(k)*W_s) + ac.P_D(t, k) - ac.C_w(k)*P_w(t) <= ...
                       ac.P_max(k)):sprintf('P_inj b%i t%i', k, t)];

                % Q_inj (2)
                C = [C, (ac.Q_min(k) ...
                    <= trace(ac.Ybar_k(k)*W_s) + ac.Q_D(t, k) <=...
                       ac.Q_max(k)):sprintf('Q_inj b%i t%i', k, t)];
            end
            % V_bus (3)
            C = [C, (ac.V_min(k)^2 ...
                <= trace(ac.M_k(k)*W_s) <= ...
                   ac.V_max(k)^2):sprintf('V_bus b%i t%i', k, t)];
            

        end

        for l = 1:ac.N_l
            % apparent line flow constraints
            if not(isa(ac.lineflows(l, W_s), 'double'))
                C = [C, (ac.lineflows(l, W_s) <= 0):sprintf('lineflows l%i t%i', l, t)];
            end
        end

        % refbus constraints
        if not(isempty(ac.refbus))
            refbus = ac.refbus + ac.N_b;    
            C = [C, (W_s(refbus, refbus) == 0):sprintf('refbus t%i', t)];
        end
    else
        % define constraints for every bus
        for k = 1:ac.N_b
            % P_inj (1)
            C = [C, ac.P_min(k) ...
                <= trace(ac.Y_P(k)*W_s) + ac.P_D(t, k) - ac.C_w(k)*P_w(t) <= ...
                   ac.P_max(k)];

            % Q_inj (2)
            C = [C, ac.Q_min(k) ...
                <= trace(ac.Y_Q(k)*W_s) + ac.Q_D(t, k) <=...
                   ac.Q_max(k)];

            % V_bus (3)
            C = [C, ac.V_min(k)^2 ...
                <= trace(ac.M_k(k, 1)*W_s) <= ...
                   ac.V_max(k)^2];
        end

        for l = 1:ac.N_l
            % apparent line flow constraints
            C = [C, ac.lineflows_complex(l, W_s) <= 0];
        end     

end