
classdef AC_agent < handle
    
    properties
        Obj;            % objective function
        ac;             % ac model
        wind;           % wind realizations
        C_0;            % set of initial deterministic constraints
        C_1_params;     % array with params of initial constraints
        A;              % array with params of active constraints
        L;              % set of constraints for the next iteration
        J;              % value of the objective function
        Jtilde;         % max of neigbouring objectives
        x;              % sdpvars of the decision variable
        x_hist;         % previous versions of x
        t;              % iteration number
        t_wind;         % wind time step
    end
    
    methods

        function ag = AC_agent(ac, wind, t_wind, i_start, i_end)
        % creates constraint set and objective function
            
            % INITIALIZE SDPVARS
            ag.x = {    sdpvar(2*ac.N_b), ...       Wf
                        sdpvar(2*ac.N_b), ...       Wmus
                        sdpvar(2*ac.N_b), ...       Wmds
                        sdpvar(2*ac.N_G, 1)}; ...   Rus and Rds
            
            % create objective function
            ag.Obj = AC_f(ag.x, ac, wind, t_wind);

            % create deterministic constraints
            ag.C_0 = AC_cons_det(ag.x, ac, wind, t_wind);
            
            % loop over scenarios to create scenario constraints
            C_ineqs = [];
            ag.C_1_params = [];
            
            for i = i_start:i_end        
                
                % inequality constraints
                C_ineq = AC_cons_scen(ag.x, ac, wind.slice(i), t_wind);
                
                % create parameters, last one is -1 for psd constraint
                Ncons = length(C_ineq);
                C_params = [ones(Ncons,1)*i [[1:Ncons-1]'; -1]];
                
                % store params
                ag.C_1_params = [ag.C_1_params; C_params];
                
                % store constraints
                C_ineqs = [C_ineqs, C_ineq];

            end
            
            ag.wind = wind;
            ag.ac = ac;
            ag.t_wind = t_wind;
        
            % optimize
            opt = sdpsettings('verbose', 0);
            diagnostics = optimize([ag.C_0, C_ineqs], ag.Obj, opt);
            assert(not(diagnostics.problem), sprintf(...
                'Problem initializing agent: %s', diagnostics.info));

            % store value of objective function
            ag.J(1) = value(ag.Obj);
            ag.Jtilde = -1e9;
            
            % store params of active constraints
            ag.A{1} = [];
            x_star = values_cell(ag.x);
            
            
            % adapt this to the new AC_active function
            for i = i_start:i_end
                
                % check if the problem is indeed feasible
                problem = AC_check(x_star, ac, wind.slice(i), t_wind);
                assert(not(problem),  ...
                     sprintf('A posteriori check failed: %i', problem));
                 
                act_params = AC_active(x_star, ac, wind.slice(i), t_wind);
                
                ag.A{1} = [ag.A{1};  [ones(length(act_params),1)*i ...
                                      act_params]];
            end
            
            % set t to 1
            ag.t = 1;
            
            % init L
            ag.L = [];
            
            % store X
            ag.x_hist{1} = x_star;

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
                x_star = values_cell(ag.x);
                for params = ag.L'
                    i = params(1);
                    j = params(2);
                    
                    act_param = AC_active(x_star, ag.ac, ag.wind.slice(i), ...
                                          ag.t_wind, j);
                    
                    % check for infeasibility
                    if isempty(act_param)
                        still_feasible = 0;
                        break;
                    end
                    
                    % store the active constraints
                    ag.A{ag.t + 1} = [ag.A{ag.t+1}; ...
                                      [ones(length(act_param))*i act_param]];
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
                               AC_cons_scen(ag.x, ag.ac, ag.wind.slice(i), ...
                                            ag.t_wind, j)];
                    end
                
                    % optimize again
                    opt = sdpsettings('verbose', 0);
                    diagnostics = optimize([ag.C_0, C_L], ag.Obj, opt);
                    assert(not(diagnostics.problem), sprintf(...
                      'Problem optimizing agent: %s', diagnostics.info));
                    
                    % update active constraints using the updated solution
                    ag.A{ag.t + 1} = [];
                    x_star = values_cell(ag.x);
                    for params = ag.L'
                        i = params(1);
                        j = params(2);
                        
                        act_param = AC_active(x_star, ag.ac, ...
                                              ag.wind.slice(i), ...
                                              ag.t_wind, j);

                        ag.A{ag.t+1} = [ag.A{ag.t+1}; ...
                                        [ones(length(act_param))*i act_param]]; 
                    end
                    
                    % update J
                    ag.J(ag.t + 1) = value(ag.Obj);
                    
                    
                end                

                % update iteration number
                ag.t = ag.t + 1;

                % reset L
                ag.L = [];
                    
                % store x
                ag.x_hist{ag.t + 1} = x_star;
            end
        end
        
    end
end