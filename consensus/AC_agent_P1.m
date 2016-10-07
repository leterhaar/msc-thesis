%% AC_AGENT_P1

classdef AC_agent_P1 < handle
    
    % all the properties containing SDP vars
    properties (Transient = true)
        cons;           % Constraints
        obj;            % Objective function
        
    end
        
    properties (Transient = false)
        ac;             % AC model
        wind;           % wind model (or batch of it)
        t;              % current timestep
        info;           % information about solving per iteration
        k;              % current iteration number
        solution;       % struct containing the solution
    end
    
    methods
       
        function ag = AC_agent_P1(ac, wind, t)
        %% Initialize
        % 
        % Arguments
        % ---------
        % ac : instance of AC_model
        % wind : instance of wind_model or struct with fields P_w, P_wf, P_m
        % t : current timestep
        
            ag.ac = ac;
            ag.wind = wind;
            ag.t = t;
            
            N = size(wind.P_w, 2);

            % Initialize decision variables
            W_f = sdpvar(2*ac.N_b);                 % forecasted voltage vectors squared
            W_s = sdpvar(2*ac.N_b, 2*ac.N_b, N);    % scenario voltage vectors squared
            R_us = sdpvar(ac.N_G, 1);               % upspinning reserve requirement
            R_ds = sdpvar(ac.N_G, 1);               % downspinning reserve requirement
            d_us = sdpvar(ac.N_G, 1);               % upspinning distribution vector
            d_ds = sdpvar(ac.N_G, 1);               % downspinnign distribution vector
            
            
            % Define objective function
            
            
            % Define constraints
            
            
            
            
            % iteration number
            ag.k = 0;
            
            % initiate result struct
            ag.solution = struct('W_f', {}, 'R_us', {}, 'R_ds', {}, 'd_us', {}, 'd_ds', {});
            ag.info = struct('exit_code', {}, 'solver_msg', {});
        end

        function Xtilde = update(ag, Z, c)
        %% Update
        % 
        % Arguments
        % ---------
        % Z : set of global variables 
        % c : consensus parameter
        
            % update iteration number
            ag.k = ag.k + 1;
            
            % Augment agentective function with consensus term
            Objective = ag.obj + 1/(2*c) * ( norm(Z{1}-ag.dec_vars{1}, 'fro')^2 + ...
                                                       norm(Z{2}-ag.dec_vars{2}, 'fro')^2 + ...
                                                       norm(Z{3}-ag.dec_vars{3})^2 + ...
                                                       norm(Z{4}-ag.dec_vars{4})^2);

            % Optimize
            optimize_settings = sdpsettings('verbose', 0);
            optim_info = optimize(ag.cons, Objective, optimize_settings);
            ag.info(ag.k).exit_code = optim_info.problem;
            ag.info(ag.k).solver_msg = optim_info.info;
                                         
            % Extract next X
            var_names = fieldnames(ag.X);
            N_G = size(ag.dec_vars{3}, 1); % R_us is N_G x 1
            twoN_b = size(ag.dec_vars{1}, 1); % W_f is 2N_b x 2N_b
            Xtilde = NaN(2*(twoN_b^2)+2*N_G, 1);
            j = 0;
            for i = 1:4
                var_value = value(ag.dec_vars{i});
                ag.solution(ag.k).(var_names{i}) = var_value;
                var_value_vector = var_value(:);
                var_n = length(var_value_vector);
                Xtilde(j+1:j+var_n) = var_value_vector;
                j = j + var_n;
                decided_vars{i} = var_value;
            end
        end
    end
end