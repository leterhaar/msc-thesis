
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

            % create deterministic constraints
            ag.C_0 = DC_f_0(ag.x, dc, wind, t_wind);
            
            % loop over scenarios to create scenario constraints
            C_ineqs = [];
            ag.C_1_params = [];
            
            for i = i_start:i_end                    
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
            optimize([ag.C_0, C_ineqs], ag.Obj, opt);

            % store value of objective function
            ag.J(1) = value(ag.Obj);
            ag.Jtilde = -1e9;
            
            % store params of active constraints
            ag.A{1} = [];
            x_star = value(ag.x);
            for i = i_start:i_end
                ag.A{1} = [ag.A{1};  ...
                       DC_f_check(x_star, i, dc, wind, t_wind)];
            end
            
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

                % add own initial constraints and last active constraints
                ag.L = unique([ag.L; ag.C_1_params; ag.A{ag.t}], 'rows');
                
                % check feasibility and activeness for all constraints
                still_feasible = 1;
                ag.A{ag.t + 1} = [];
                x_star = value(ag.x);
                for params = ag.L'
                    i = params(1);
                    j = params(2);
                    
                    [params_act, residuals] = DC_f_check(x_star, i, ag.dc,...
                                                ag.wind, ag.t_wind, j);
                    
                    % check for infeasibility
                    if residuals < -1e-6 || isnan(residuals)
                        still_feasible = 0;
                        break;
                    end
                    
                    % store the active constraints
                    ag.A{ag.t + 1} = [ag.A{ag.t+1}; params_act];
                end
                
                % see if the solution is still feasible
                if still_feasible
                    
                    % keep the objective value the same
                    ag.J(ag.t + 1) = ag.J(ag.t);
                    
                else
                    
                    % build constraints from L
                    C_L = [];
                    for params = ag.L'

                        % extract params
                        i = params(1);
                        j = params(2);

                        % add constraints to set C_L
                        C_L = [C_L, ...
                               DC_f_ineq(ag.x, i, ag.dc, ag.wind, ag.t_wind, j)];
                    end
                
                    % optimize again
                    opt = sdpsettings('verbose', 0);
                    optimize([ag.C_0, C_L], ag.Obj, opt);
                    
                    % update active constraints using the updated solution
                    ag.A{ag.t + 1} = [];
                    x_star = value(ag.x);
                    for params = ag.L'
                        i = params(1);
                        j = params(2);
                        ag.A{ag.t+1} = [ag.A{ag.t+1}; ...
                                        DC_f_check(x_star, i, ag.dc, ...
                                        ag.wind, ag.t_wind, j)];
                    end
                    
                    % update J
                    ag.J(ag.t + 1) = value(ag.Obj);
                    
                end                

                % update iteration number
                ag.t = ag.t + 1;

                % reset L
                ag.L = [];
            end
        end
        
    end
end