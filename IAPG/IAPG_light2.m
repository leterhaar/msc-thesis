%% [xstar, iterations] = IAPG(x, f, gradient_f, constraints, ...)
% runs the incremtal aggregated proximal gradient algorithm and returns 
% solution and iterations
% 
% solves min_ f(x) s.t. x \in X_i for i = 1,...,m
%
% PARAMETERS
% * x               : an sdpvar
% * fs              : a cell of m function handles to the smooth part of 
%                     the objective function f_i(x)
% * grad_fs         : a cell of m function handles to the gradient of the smooth part
%                     of the objective function gradient_f(x)
% * constraints     : a cell of constraint sets
% * options         : set of options, optional, as follows
%   - x0            : value for x0. if empty, a feasible x0 will be found
%   - verbose       : 0 (default) or 1 (show progress)
%   - max_its       : maximum number of iterations (default 500)
%   - opt_settings  : optimization settings (default verbose = 0)
%   - default_constraint : default constraint for every optimization step
%   - alpha         : function handle for alpha. default alpha = 1/(k+1)


function [xstar, its] = IAPG_light2(x_sdp, fs, grad_fs, constraints, varargin)

    %% Load options and check validity of inputs
    
    % check types of input
    assert(isa(x_sdp, 'sdpvar'), 'x should be sdpvar');
    assert(isa(fs, 'cell') && isa(fs{1}, 'function_handle'), ...
                                'fs should be cell with function handles');
    assert(isa(grad_fs, 'cell') && isa(grad_fs{1}, 'function_handle'), ...
                           'grad_fs should be cell with function handles');
    assert(isa(constraints, 'cell') && (isa(constraints{1}, 'constraint') || ...
      isa(constraints{1}, 'lmi')), 'constraints should be cell with lmis');
    
    % define default options
    options = struct('verbose', 0, 'max_its', 500, 'opt_settings', ...
                     sdpsettings('verbose', 0), 'x0', [], ...
                     'default_constraint', [], 'alpha', @(k) 1/(k+1));
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
    
    % store dimensions
    d = size(x_sdp);
    m = length(constraints);
    
    %% Initialize algorithm
    
    % find feasible x0
    if isempty(options.x0)
        info = optimize([options.default_constraint, ...
                            [constraints{:}]], [], options.opt_settings);
        assert(not(info.problem), 'Problem optimizing: %s', info.info);
        x0 = value(x_sdp);
    else
        assert(all(size(options.x0) == d), 'x0 has the wrong size');
        x0 = options.x0;
    end
    
    its = struct(   'x', [], ...        value for x
                    'f', [], ...        objective function value
                    'subgrad', [], ...  subgradient for projection at x
                    'time', []...       time per iterations    
                    );
        
    its(1).x = x0;
                
    %% cache solvers
    
    if verbose
        p = progress('Preparing IAPG-light2', m);
    end
    
    if isempty(options.opt_settings.solver)
        warning('Calling optimizer with no solver specified!');
    end
    
    % create sdpvar for z, past subgrads, ak
    z_sdp = sdpvar(d(1), d(2), 'full');
    a_sdp = sdpvar(1);
    
    % define objective with vectorized versions (so matrix also works)
    Obj = 1/(2*a_sdp) * norm(x_sdp(:) - z_sdp(:), 2)^2;
    
    proximals = cell(m,1);
    for i = 1:m
        proximals{i} = optimizer([options.default_constraint; ...
                                 [constraints{i}]], Obj,  ...
                                  options.opt_settings, ...
                                  {a_sdp, z_sdp}, ...
                                  x_sdp);
                              
        if verbose
            p.ping();
        end
    end
    
    %% prepare first m iterations

    for k = 1:2*m
        % loop over all agents
        i = rem(k-1, m) + 1;
        ak = options.alpha(k);
        
        % gradient step (no aggregation)
        zk = its(k).x - ak * grad_fs{i}(its(k).x);
        
        % proximal step
        [next_x, problem, msg] = proximals{i}(ak, zk);
        assert(not(problem), sprintf('Problem optimizing: %s', msg{:}));
        
        its(k+1).x = next_x;
        its(k+1).f = fs{i}(next_x);
        its(k+1).grad = grad_fs{i}(next_x);
        its(k+1).subgrad = 1/ak * (zk - next_x);
        
    end
    
    % store the most recent iterations for every agent
    li_k = m+1:2*m;

    %% start main iterations

    if verbose
        p = progress('Iterating IAPG light2', options.max_its-2*m);
    end
    
    all_is = 1:m;
    
    for k = 2*m:options.max_its
        tic
        
        xk = its(k).x;
        ak = options.alpha(k);
        
        % select a random agent
        i = randi(m);           % random
%         i = rem(k-1, m) + 1;  % cycling
        all_except_i = all_is(all_is ~= i);

        % sum past (sub)gradients
        past_gradients = zeros_like(x0);
        
        for l_i = li_k(all_except_i)
            
            % sum gradients for all agents except current agent
            past_gradients = past_gradients + its(l_i).grad;
            
        end
        
        zk = xk - ak * (grad_fs{i}(xk) + past_gradients);
        
        [next_x, problem, msg] = proximals{i}(ak, zk);
        assert(not(problem), sprintf('Problem optimizing: %s', msg{:}));

        % store x(k+1) and its gradient and objective function
        its(k+1).x = next_x;
        its(k+1).f = fs{i}(next_x);
        its(k+1).time = toc;
        its(k).grad = grad_fs{i}(xk);

        li_k(i) = k;
        
        if verbose
            p.ping(); 
        end
    
    end
    
    % return the final value and iterations
    xstar = next_x;
%     its = its(2*m:end); % remove first 2m, only used for initialization
end