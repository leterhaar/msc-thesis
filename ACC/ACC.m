%% [xstar, its] = ACC(x_sdp, delta_sdp, deltas, f, constraints, ...)
% Executes the Active Constraint Consensus algorithm to solve the
% optimization problem: 
%  min_x f(x)
%  s.t. constraints(x, delta_i) for i = 1, ..., N
% 
% Uses the optimizer functionality from YALMIP
%
% PARAMETERS
% ==========
% x_sdp         : optimization variable (sdpvar)
% delta_sdp     : parameter variable (sdpvar)
% f             : function handle for the objective function
% constraints   : LMI object with all the constraints as a function of
%                 x_sdp and delta_sdp
% deltas        : a N x .. matrix with realizations of delta on the rows
% options       : key value pairs, optional, of following format:
%  - verbose    : flag to show progress, 1 = show (default), 0 = hide
%  - opt_settings : sdpsettings for optimization (default verbose=0)
%  - default_constriant : default (deterministic constraint)
%  - diameter   : diameter of connectivity graph (default 3)
%  - n_agents   : number of agents (default ceil(N/10))
%  - debug      : enter debugging inside function on error (default 0)
%  - x0         : initial value for x (if empty, zeros)
%  - max_its    : maximum no of iterations
%
% RETURNS
% =======
% xstar         : optimal value for x after convergence
% agents        : structure with iterations

function [xstar, agents] = ACC(x_sdp, delta_sdp, deltas, f, constraints, varargin)
    %% check validity of input and load options
    % check types of input
    assert(isa(x_sdp, 'sdpvar'), 'x_sdp should be sdpvar');
    assert(isa(delta_sdp, 'sdpvar'), 'delta_sdp should be sdpvar');
    assert(isa(f, 'function_handle'), ...
                                'f should be function handle');
    assert(isa(constraints, 'constraint') || ...
      isa(constraints, 'lmi'), 'constraints should be lmis');
  
    % store dimensions
    d = size(x_sdp);
    [N, d_delta] = size(deltas);
    
    % check dimensions
    assert(length(delta_sdp) == d_delta);
    
    % define default options
    options = struct('verbose', 0, ...
                     'opt_settings', sdpsettings('verbose', 0), ...
                     'x0', [], ...
                     'default_constraint', [], ...
                     'diameter', 3, ...
                     'debug', 0, ...
                     'n_agents', [],...
                     'max_its', 15);
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
        options.x0 = zeros_like(x_sdp);
    end
    
    
    
    try
    
        %% create connectivity graph
        if not(exist('random_graph', 'file'))
            % add to path
            addpath('../misc');
        end

        connectivity_graph = random_graph(m,options.diameter, 'rand');
        
        %% initialize agents
        agents = struct('initial_deltas', [], ...
                        'iterations', []);
        
        % divide the initial deltas among the agents
        Nm = ceil(N/m); % no of constraints per agent
        for i = 1:m
            agents(i).initial_deltas = deltas(((i-1)*Nm)+1:min(i*Nm,N), :);
        end
        
        [agents.iterations] = deal(struct(  'J', nan, ...
                                            'active_deltas', [],...
                                            'x', options.x0, ...
                                            'time', nan));
                  
        % precompile solver
        the_solver = optimizer([options.default_constraint, constraints], ...
                                f(x_sdp), options.opt_settings, delta_sdp, x_sdp);
                        
        %% start main iterations
        k = 1;
        ngc = ones(m, 1);
        
        while all(ngc < 2*options.diameter+1) && k < options.max_its
            
            if verbose
                prg = progress(sprintf('Iteration %i',k), m);
            end

            % loop over agents
            for i = 1:m
                
                tic
                
                % build constraint set from C_i and A_i ...
                L = [agents(i).initial_deltas; ...
                     agents(i).iterations(k).active_deltas];
                 
                % ... and A_j for all incoming agents j
                for j = find(connectivity_graph(:, i))';
                    L = [L; agents(j).iterations(k).active_deltas];
                end

                % filter out double deltas
                L = unique(L, 'rows');
                
                % build solver
                merged = [];
                feasible_for_all = 1;
                for j = 1:size(L,1)
                    
                    % check feasibility of new set of constraints
                    assign(x_sdp, agents(i).iterations(k).x);
                    assign(delta_sdp, L(j, :));
                    residuals = check(constraints);
                    if any(residuals < -1e-6);
                        feasible_for_all = 0;
                    end
                    
                    % merge the solvers with filled deltas together
                    merged = [merged; the_solver(L(j,:), 'nosolve')];
                end
                
                if not(feasible_for_all) || k == 1

                    % solve
                    [next_x, problem, msg] = merged([]);
                    assert(not(problem), 'Problem optimizing: %s', msg{:});
                    
                    agents(i).iterations(k+1).x = next_x;
                    agents(i).iterations(k+1).J = f(next_x);
                    
                    
                    % store active constraints
                    agents(i).iterations(k+1).active_deltas = [];
                    for j = 1:size(L,1)
                         % check feasibility of new set of constraints
                        assign(x_sdp, next_x);
                        assign(delta_sdp, L(j, :));
                        residuals = check(constraints);
                        if any(residuals < 1e-6 & residuals > -1e-6);
                            agents(i).iterations(k+1).active_deltas = [
                                agents(i).iterations(k+1).active_deltas;
                                L(j, :)];
                        end
                    end
                else
                    agents(i).iterations(k+1).x = ...
                                                agents(i).iterations(k).x;
                    agents(i).iterations(k+1).J = ...
                                                agents(i).iterations(k).J;

                end
                
                % check if J(t+1) is J(t)
                if all_close(agents(i).iterations(k).J, agents(i).iterations(k+1).J, 1e-6)
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
        end
        
        % check if everything went well
        for i = 1:m
            for j = i+1:m
                assert(all_close(agents(i).iterations(k).x, ...
                                 agents(j).iterations(k).x, 1e-3), ...
                                 'agents not close'); 
            end
        end
        
        xstar = agents(1).iterations(k).x;
            

        %% return output
    catch e
        if options.debug==1
            fprintf('\n%s', e.getReport());
            keyboard;
        else
            rethrow(e);
        end
    end
end