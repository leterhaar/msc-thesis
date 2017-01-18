%% [xstar, its] = ACCA(x_sdp, deltas, objective_fcn, constraints_fcn, ...)
% Executes the Active Constraint Consensus Agreement algorithm to solve the
% optimization problem: 
%  min_x f(x)
%  s.t. constraints(x, delta_i) for i = 1, ..., N
%
% PARAMETERS
% ==========
% x_sdp         : cell with sdpvars
% deltas        : a N x .. matrix with realizations of delta on the rows
% objective_fcn : function handle for the objective function
% constraints_fcn : a function that returns an LMI object 
% options       : key value pairs, optional, of following format:
%  - verbose    : flag to show progress, 1 = show (default), 0 = hide
%  - opt_settings : sdpsettings for optimization (default verbose=0)
%  - default_constriant : default (deterministic constraint)
%  - diameter   : diameter of connectivity graph (default 3)
%  - n_agents   : number of agents (default ceil(N/10))
%  - debug      : enter debugging inside function on error (default 0)
%  - x0         : initial value for x (if empty, zeros)
%  - max_its    : maximum no of iterations (default 100)
%  - residuals  : function handle to evaluate residuals h(x,delta) >= 0
%                 optional, when empty the check function of yalmip will be
%                 used
%  - use_selector : boolean to indicate that the constraint function
%                   accepts a third argument h(x, delta, j) >= 0 for
%                   selection. 
%  - connectivity : adjacency matrix of the connectivity graph
%  - stepsize   : stepsize function handle of the form @(k) ...,          
%                 default 1/(k+1)
%  - tolerance  : tolerance for determinining activeness (default 1e-6)
%
% RETURNS
% =======
% xstar         : optimal value for x after convergence
% agents        : structure with iterations

function [xstar, agents] = ACCA_fcn_cell(x_cell, deltas, f, cons_fcn, varargin)
    %% check validity of input and load options
    % check types of input
    assert(isa(x_cell, 'cell') && isa(x_cell{1}, 'sdpvar'), 'x_sdp should be sdpvar');
    assert(isa(f, 'function_handle'), 'f should be function handle');
    assert(isa(cons_fcn, 'function_handle'), ...
                                       'constraints should be a function');
    
    N = size(deltas, 1);
    N_cell = length(x_cell); % number of elements in the cell

    % define default options
    options = struct('verbose', 0, ...
                     'opt_settings', sdpsettings('verbose', 0), ...
                     'x0', [], ...
                     'default_constraint', [], ...
                     'diameter', 3, ...
                     'debug', 0, ...
                     'n_agents', [],...
                     'max_its', 100, ...
                     'residuals', [],...
                     'use_selector', false,...
                     'connectivity', [],...
                     'stepsize', @(k)1/(k+1),...
                     'tolerance', 1e-6);
    def_fields = fieldnames(options);
                 
    % load options from varargin
    nvarargin = length(varargin);
    assert(rem(nvarargin,2) == 0, 'Please provide key-value pairs');
    fields = varargin(1:2:end);
    values = varargin(2:2:end);
    
    % override options 
    for i = 1:length(fields)
        if find(strcmp(def_fields, fields{i}))
            options.(fields{i}) = values{i};
        else
            warning('Field "%s" is unknown. Typo?', fields{i});
        end
    end
                 
    % set verbose flag
    if options.verbose == 1
        verbose = true;
    else
        verbose = false;
    end
    
    if options.debug == 1
        debug = true;
    else
        debug = false;
    end
    
    % set number of agents
    if isempty(options.n_agents)
        m = ceil(N/10);
    else
        assert(isnumeric(options.n_agents), 'Should be numeric');
        m = options.n_agents;
    end
    assert(N >= m, ...
            'Number of agents can not be larger than number of scenarios');
    
    
    if isempty(options.opt_settings.solver)
        warning('Calling optimizer with no solver specified!');
    end
    
    % add helper path
    if not(exist('zeros_like', 'file'))
        addpath('../misc/');
    end
    
    if isempty(options.x0)
        options.x0 = cell(N_cell, 1);
        for i = 1:N_cell
            options.x0{i} = zeros_like(x_cell{i});
        end
    end
    
    % find connectivity graph
    if isempty(options.connectivity)
        if not(exist('digraph', 'file'))
            warning('Digraph function not found, using G from workspace');
            connectivity_graph = evalin('base', 'G');
        else 
            connectivity_graph = random_graph(m,options.diameter, 'rand');
        end
    else
        connectivity_graph = options.connectivity;
    end
    
    % use try block to enter debug mode inside function when encountering
    % error
    try      
        
        %% initialize agents
        agents = struct('initial_deltas', [], ...
                        'iterations', []);
        
        % find the number of constraints that results from the cons_fcn
        if options.use_selector
            Ncons = length(cons_fcn(x_cell, deltas(1, :), 0));
        else
            Ncons = length(cons_fcn(x_cell, deltas(1, :)));
        end

        % divide the initial deltas among the agents
        Nm = ceil(N/m);             % no of deltas per agent
        constraint_ids = [1:Ncons]';  % list with constraint identifiers
        for i = 1:m
            % store the deltas in the agent; the first entry of the row
            % will be the identifier of the constraint, the rest of the row
            % will be the actual data corresponding to that delta
            delta_slice = deltas(((i-1)*Nm)+1:min(i*Nm,N), :);
            identifiers = repmat(constraint_ids, size(delta_slice, 1), 1);
            agents(i).initial_deltas = [identifiers, ...
                                        kron(delta_slice, ones(Ncons, 1))];
                
        end
        
        [agents.iterations] = deal(struct(  'J', f(options.x0), ...
                                            'active_deltas', [],...
                                            'x', {options.x0}, ...
                                            'time', nan,...
                                            'info', struct('num_cons', 0, ...
                                                           'optimized', 1)));
                       
        %% start main iterations
        k = 1;
        ngc = ones(m, 1);
        loop_active = 1;
        
        while loop_active
            
            if verbose
                prg = progress(sprintf('Iteration %i',k), m);
            end

            % loop over agents
            for i = 1:m
                
                tic
                %% build C_i and z
                
                % build constraint set from C_i and A_i ...
                L = [agents(i).initial_deltas; ...
                     agents(i).iterations(k).active_deltas];
                 
                % build consensus variable z from average of own solution +
                % incoming constraints using constant weights a_*^i =
                % 1/(N+1)
                N_incoming = sum(connectivity_graph(:,i)) + 1;
                z = cell_divide(agents(i).iterations(k).x, N_incoming);
                
                % loop over neighbouring agents
                for j = find(connectivity_graph(:, i))';
                    % add A_j for all incoming agents j to constraint set
                    L = [L; agents(j).iterations(k).active_deltas];
                    
                    % sum up incoming xs devided by number of connections
                    z = cell_add(z,cell_divide(agents(j).iterations(k).x,N_incoming));
                end

                % filter out double deltas
                L = unique(L, 'rows');
                
                %% check feasibility of new set of constraints
                
                feasible_for_all = 1;
                for j = 1:size(L,1)
                    
                    if isempty(options.residuals) % use YALMIP check
                        cons_delta = cons_fcn(x_cell, L(j, 2:end));
                        assign_cell(x_cell, agents(i).iterations(k).x);
                        residual = check(cons_delta(L(j, 1)));
                        
                    elseif options.use_selector % use h(x, delta, j) >= 0
                        residual = options.residuals(...
                                                agents(i).iterations(k).x, ...
                                                L(j, 2:end), L(j, 1));
                                            
                    else % use residual function h(x, delta) >= 0
                        % get all residuals
                        residuals = options.residuals(...
                                                agents(i).iterations(k).x, ...
                                                L(j, 2:end));
                        % filter out the residual of interest
                        residual = residuals(L(j,1));
                    end
                    
                    if residual < -options.tolerance
                        feasible_for_all = 0;
                        break
                    end
                end
                
                %% update x
                
                % optimize if new z is not feasible
                if not(feasible_for_all) || k == 1
                    
                    N_cons_used = 0;
                    % build solver
                    C_i = [];
                    C_all = options.default_constraint; 
                    for j = 1:size(L,1)
                        
                        % check feasibility for z instead of x
                        if isempty(options.residuals) % use YALMIP check
                            cons_delta = cons_fcn(x_cell, L(j, 2:end));
                            assign_cell(x_cell, z);
                            residual = check(cons_delta(L(j, 1)));

                        elseif options.use_selector % use h(x, delta, j) >= 0
                            residual = options.residuals(z, L(j, 2:end),...
                                                         L(j, 1));

                        else % use residual function h(x, delta) >= 0
                            % get all residuals
                            residuals = options.residuals(z, L(j, 2:end));
                            % filter out the residual of interest
                            residual = residuals(L(j,1));
                        end

                        % only use the constraints from infeasible or
                        % strict deltas
                        if residual < options.tolerance % || true % REMOVE '|| true' FOR FULL ACCA
                            C_i = [C_i; L(j, :)];
                            
                             % use previously defined cons_delta
                            if isempty(options.residuals)
                                C_all = [C_all, cons_delta(L(j,1))];
                                
                            % use cons_fcn with selector
                            elseif options.use_selector
                                C_all = [C_all, cons_fcn(x_cell, ...
                                                             L(j, 2:end),...
                                                             L(j, 1))];
                            % use cons_fcn without selector
                            else
                                cons_delta = cons_fcn(x_cell, L(j, 2:end));
                                C_all = [C_all, cons_delta(L(j, 1))];
                                
                            end
                                           
                            N_cons_used = N_cons_used + 1;
                        end
        
                    end
                    
                    if N_cons_used == 0
                        warning('Not using any scenario constraints...');
                        % we should not get here....
                    end

                    % define objective and solve
                    alpha = options.stepsize(k);
                    Obj_consensus = f(x_cell);
                    % UNCOMMENT BELOW FOR FULL ACCA
                    for l = 1:N_cell
                        Obj_consensus = Obj_consensus + 1/(2*alpha) * ...
                                            norm(x_cell{l}(:) - z{l}(:))^2;
                    end
                    % UNCOMMENT ABOVE FOR FULL ACCA
                    status = optimize(C_all, Obj_consensus, ...
                                                     options.opt_settings);
                    assert(not(status.problem), status.info);
                    next_x = values_cell(x_cell);
                    
                    agents(i).iterations(k+1).x = next_x;
                    agents(i).iterations(k+1).J = f(next_x);
                    
                    % store status
                    if debug
                        agents(i).iterations(k+1).info.optimized = 1;
                        agents(i).iterations(k+1).info.num_cons = N_cons_used;
                    end
                    
                else
                    agents(i).iterations(k+1).x = ...
                                                agents(i).iterations(k).x;
                    agents(i).iterations(k+1).J = ...
                                                agents(i).iterations(k).J;
                    
                    % store status
                    if debug
                        agents(i).iterations(k+1).info.optimized = 0;
                        agents(i).iterations(k+1).info.num_cons = size(L,1);
                    end
                end
                
                % store active constraints
                agents(i).iterations(k+1).active_deltas = [];
                assign_cell(x_cell, agents(i).iterations(k+1).x);

                for j = 1:size(L,1)
                    
                    % check feasibility of of the new solution
                    if isempty(options.residuals) % use YALMIP check
                        cons_delta = cons_fcn(x_cell, L(j, 2:end));
                        residual = check(cons_delta(L(j, 1)));
                        
                    elseif options.use_selector % use h(x, delta, j) >= 0
                        residual = options.residuals(...
                                                agents(i).iterations(k+1).x, ...
                                                L(j, 2:end), L(j, 1));
                                            
                    else % use residual function h(x, delta) >= 0
                        % get all residuals corresponding to the delta
                        residuals = options.residuals(...
                                                agents(i).iterations(k+1).x, ...
                                                L(j, 2:end));
                        % filter out the residual of interest
                        residual = residuals(L(j,1));
                    end

                    if residual < options.tolerance && residual > -options.tolerance
                        agents(i).iterations(k+1).active_deltas = [
                            agents(i).iterations(k+1).active_deltas;
                            L(j, :)];
                    end
                end
                
                % check if J(t+1) is J(t)
                if all_close(agents(i).iterations(k).J, agents(i).iterations(k+1).J, options.tolerance)
                    ngc(i) = ngc(i) + 1;
                else
                    ngc(i) = 1;
                end
                
                if verbose
                    prg.ping();
                end
                
                agents(i).iterations(k+1).time = toc;
            end

            % update iteration number     
            k = k + 1;
            
            if all(ngc >= 2*options.diameter+1)
                loop_active = 0;
            elseif k >= options.max_its
                warning(['Maximum iterations (%i) is reached, '...
                           'might not have convergence'], options.max_its);
                loop_active = 0;
            end
        end
        
        % check if everything went well
        for i = 1:m
            for j = i+1:m
                if not(all_close(agents(i).iterations(k).x, ...
                                 agents(j).iterations(k).x, 1e-3))
                    warning('Agents not close, xstar may not be optimal'); 
                end
            end
        end
        
        xstar = agents(1).iterations(k).x;
            

        %% return output
    catch e
        if debug
            fprintf('\n%s', e.getReport());
            keyboard;
            xstar = agents(1).iterations.x;
        else
            rethrow(e);
        end
    end
end