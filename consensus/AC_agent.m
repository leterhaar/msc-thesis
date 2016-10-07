%% DC Agent
% An agent in the distrituted consensus algorithm

classdef AC_agent < handle
   
    properties (Transient = false)
        X;              % struct with agent local variables
        problem;        % problem formulation with local constraints
        info;           % struct with history of past optimizations
        k;              % iteration number
    end
    
    properties (Transient = true)
        dec_vars;       % decision variables (sdpvars)
    end

    methods
        
        function agent = AC_agent(ac, wind, t)
        %% Initialize
        % 
        % Arguments
        % ---------
        % dc : instance of DC_model
        % wind : instance of wind_model or struct with field P_w, P_wf, P_m

            % Initialize decision variables
            W_f = sdpvar(2*ac.N_b);                % forecasted voltage vectors squared
            W_m = sdpvar(2*ac.N_b);                % mismatch voltage vectors squared
            R_us = sdpvar(ac.N_G, 1);      % upspinning reserve requirement
            R_ds = sdpvar(ac.N_G, 1);      % downspinning reserve requirement
            agent.dec_vars = {W_f, W_m, R_us, R_ds};
            
            % Prepare agentective function and constraints
            agent.problem = formulate('P3*');
            agent.problem.prepare(ac, wind, t, agent.dec_vars);
            
            % iteration number
            agent.k = 0;
            
            % initiate result struct
            agent.X = struct('W_f', {}, 'W_m', {}, 'R_us', {}, 'R_ds', {});
            agent.info = struct('exit_code', {}, 'solver_msg', {});
        end

        function Xtilde = update(agent, Z, c)
        %% Update
        % 
        % Arguments
        % ---------
        % Z : set of global variables 
        % c : consensus parameter
        
            % update iteration number
            agent.k = agent.k + 1;
            
            % Augment agentective function with consensus term
            Objective = agent.problem.Objective + 1/(2*c) * ( norm(Z{1}-agent.dec_vars{1}, 'fro')^2 + ...
                                                       norm(Z{2}-agent.dec_vars{2}, 'fro')^2 + ...
                                                       norm(Z{3}-agent.dec_vars{3})^2 + ...
                                                       norm(Z{4}-agent.dec_vars{4})^2);

            % Optimize
            optimize_settings = sdpsettings('verbose', 0);
            optim_info = optimize(agent.problem.Constraints, ...
                                             Objective, optimize_settings);
            agent.info(agent.k).exit_code = optim_info.problem;
            agent.info(agent.k).solver_msg = optim_info.info;
                                         
            % Extract next X
            var_names = fieldnames(agent.X);
            N_G = size(agent.dec_vars{3}, 1); % R_us is N_G x 1
            twoN_b = size(agent.dec_vars{1}, 1); % W_f is 2N_b x 2N_b
            Xtilde = NaN(2*(twoN_b^2)+2*N_G, 1);
            j = 0;
            for i = 1:4
                var_value = value(agent.dec_vars{i});
                agent.X(agent.k).(var_names{i}) = var_value;
                var_value_vector = var_value(:);
                var_n = length(var_value_vector);
                Xtilde(j+1:j+var_n) = var_value_vector;
                j = j + var_n;
                decided_vars{i} = var_value;
            end

        end
    end
end