%% DC Agent
% An agent in the distrituted consensus algorithm

classdef DC_agent < handle
   
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
        
        function agent = DC_agent(dc, wind, t)
        %% Initialize
        % 
        % Arguments
        % ---------
        % dc : instance of DC_model
        % wind : instance of wind_model or struct with field P_w, P_wf, P_m

            % Initialize decision variables
            agent.dec_vars = sdpvar(5*dc.N_G, 1, 'full');   
            
            % split into 5 separate values for problem formulation
            problem_vars = cell(5,1);
            for i = 1:5
                var_index = (i-1)*dc.N_G + 1:i*dc.N_G;
                problem_vars{i} = agent.dec_vars(var_index);
            end
            
            % Prepare agentective function and constraints
            agent.problem = formulate('DC');
            agent.problem.prepare(dc, wind, t, problem_vars);
            
            % iteration number
            agent.k = 0;
            
            % initiate result struct
            agent.X = struct('P_G', {}, 'R_us', {}, 'R_ds', {}, 'd_us', {}, 'd_ds', {});
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
            Objective = agent.problem.Objective + ...
                        1/(2*c) * (norm(agent.dec_vars - Z)^2);

            % Optimize
            optimize_settings = sdpsettings('verbose', 0);
            optim_info = optimize(agent.problem.Constraints, ...
                                             Objective, optimize_settings);
            agent.info(agent.k).exit_code = optim_info.problem;
            agent.info(agent.k).solver_msg = optim_info.info;
                                         
            % Extract next X
            var_names = fieldnames(agent.X);
            N_G = size(agent.dec_vars, 1)/5; % var length is 5 * N_G
            Xtilde = NaN(5*N_G,1);
            for i = 1:5
                var_index = (i-1)*N_G + 1:i*N_G;
                agent.X(agent.k).(var_names{i}) = value(agent.dec_vars(var_index));
                Xtilde(var_index) = value(agent.dec_vars(var_index));
            end

        end
    end
end