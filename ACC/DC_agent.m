
classdef DC_agent < handle
    
    properties
        Obj;            % objective function
        dc;             % dc model
        wind;           % wind realizations
        C_0;            % set of initial deterministic constraints
        C_1_params;     % array with params of initial constraints
        A;              % array with params of active constraints
        L;              % set of constraints for the next iteration
        J;              % value of the objective function
        Jtilde;         % max of neigbouring objectives
        x;              % sdpbars of the decision variable
        t;              % iteration number
        t_wind;         % wind time step
    end
    
    methods
        
        function ag = DC_agent(dc, wind, t_wind, i_start, i_end)
        % creates constraint set and objective function
            
            % INITIALIZE SDPVARS
            ag.x = sdpvar(5*dc.N_G, 1, 'full'); % P_G Rus Rds dus dds
            
            % create objective function
            ag.Obj = DC_f_obj(ag.x, dc, wind, t_wind);

            % create deterministic constraint
            ag.C_0 = DC_f_0(ag.x, dc, wind, t_wind);
            
            % loop over scenarios to create scenario constraints
            C_ineqs = [];
            C_eqs = [];
            ag.C_1_params = [];
            
            for i = i_start:i_end
                % equality constraints
                C_eqs = [C_eqs, DC_f_eq(ag.x, i, dc, wind, t_wind)];
                    
                % inequality constraints
                [C_ineq, C_params] = DC_f_ineq(ag.x, i, dc, wind, t_wind);
                C_ineqs = [C_ineqs, C_ineq];
                
                % store params to inequality constraints
                ag.C_1_params = [ag.C_1_params; C_params];

            end
            
            ag.wind = wind;
            ag.dc = dc;
            ag.t_wind = t_wind;
        
            % optimize
            opt = sdpsettings('verbose', 0);
            optimize([ag.C_0, C_eqs, C_ineqs], ag.Obj, opt);

            % store value of objective function
            ag.J(1) = value(ag.Obj);
            ag.Jtilde = -1e9;
            
            % store active constraints
            ag.A{1} = ag.C_1_params(Ac(C_ineqs), :);

            % set t to 1
            ag.t = 1;
            
            % init L
            ag.L = [];

        end
        
        function build(ag, A_incoming, J_incoming)
        % builds L(t+1) and tilde J(t+1) agent by agent
            
            % add incoming A to the L(t+1)
            ag.L = [ag.L; A_incoming];
            
            % take maximum of current J(t+1) and incoming J
            ag.Jtilde = max(ag.Jtilde, J_incoming);
        
        end
        
        function update(ag)
        % optimize, update active constraint and objective function 
        
            % check feasibility neighbours
            if isinf(ag.Jtilde)
                ag.A{ag.t + 1} = [];
                ag.J{ag.t + 1} = Inf;
                ag.L = [];
            else

                % add own last active constraints and initial constraints
                ag.L = unique([ag.L; ag.C_1_params; ag.A{ag.t}], 'rows');

                % build constraints from L
                C_ineqs = [];
                for params = ag.L'
                    
                    % extract params
                    i = params(1);
                    j = params(2);
                    
                    % add constraints to set C_L
                    C_ineqs = [C_ineqs, ...
                           DC_f_ineq(ag.x, i, ag.dc, ag.wind, ag.t_wind, j)];
                end
                
                % build equality constraints
                C_eqs = [];
                for i = unique(ag.L(:,1))'
                    C_eqs = [C_eqs, ...
                                DC_f_eq(ag.x, i, ag.dc, ag.wind, ag.t_wind)];
                end

                % check if current value for x is infeasible for the new 
                % constraints
                C_L = [C_eqs, C_ineqs];
                if any(check(C_L) < -1e-6) || any(isnan(check(C_L)))

                    % if infeasible with new constraints, optimize again
                    opt = sdpsettings('verbose', 0);
                    optimize([ag.C_0, C_L], ag.Obj, opt);
                end
                    
                % update A
                ag.A{ag.t + 1} = ag.L(Ac(C_ineqs), :);

                % update J
                ag.J(ag.t + 1) = value(ag.Obj);

                % update iteration number
                ag.t = ag.t + 1;

                % reset L
                ag.L = [];
            end
        end
        
    end
end