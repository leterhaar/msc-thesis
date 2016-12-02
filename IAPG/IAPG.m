%% [xstar, iterations] = IAPG(x, f, gradient_f, constraints, ...)
% runs the incremtal aggregated proximal gradient algorithm and returns 
% solution and iterations
% 
% solves min_ f(x) s.t. x \in X_i for i = 1,...,m
%
% PARAMETERS
% * x               : an sdpvar
% * f               : a function handle to the smooth part of the objective
%                     function f(x)
% * gradient_f      : a function handle to the gradient of the smooth part
%                     of the objective function gradient_f(x)
% * constraints     : a cell of constraint sets
% * options         : set of options, optional, as follows
%   - x0            : value for x0. if empty, a feasible x0 will be found
%   - verbose       : 0 (default) or 1 (show progress)
%   - max_its       : maximum number of iterations (default 500)
%   - opt_settings  : optimization settings (default verbose = 0)
%   - default_constraint : default constraint for every optimization step

function [xstar, its] = IAPG(x_sdp, f, gradient_f, constraints, varargin)

    %% Load options and check validity of inputs
    
    % check types of input
    assert(isa(x_sdp, 'sdpvar'), 'x should be sdpvar');
    assert(isa(f, 'function_handle'), 'f should be function handle');
    assert(isa(gradient_f, 'function_handle'), ...
                                   'gradient_f should be function handle');
    assert(isa(constraints, 'cell') && (isa(constraints{1}, 'constraint') || ...
           isa(constraints{1}, 'lmi')), 'constraints should be constraints');
    
    % define default options
    options = struct('verbose', 0, 'max_its', 500, 'opt_settings', ...
                     sdpsettings('verbose', 0), 'x0', [], ...
                     'default_constraint', []);
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
        x0 = value(x);
    else
        assert(all(size(options.x0) == d), 'x0 has the wrong size');
        x0 = options.x0;
    end
    
    % make sdpvar
    % x_sdp = sdpvar(d(1), d(2), options.sdpvar_type);

    its = struct(   'x', [], ...        value for x
                    'f', [], ...        objective function value
                    'grad', [], ...     gradient at x
                    'subgrad', [], ...  subgradient for projection at x
                    'time', []...       time per iterations    
                    );
                
    %% prepare first m iterations
    [its(1:m).x] = deal(x0);
    [its(1:m).f] = deal(f(x0));
    [its(1:m).grad] = deal(gradient_f(x0));
    [its(1:m).subgrad] = deal(zeros_like(x0));
    [its(1:m).time] = deal(nan);

    %% start main iterations

    if verbose
        p = progress('Iterating IAPG', options.max_its - m);
    end
    
    for k = m:options.max_its
        tic
        % select the next agent
        xk = its(k).x;
        ak = 0.005;
        % select a random agent
        i = randi(m);

        % sum past (sub)gradients
        past_gradients = zeros_like(x0);
        past_subgrads = zeros_like(x0);
        
        for l = 0:m-1
            past_gradients = past_gradients + its(k-l).grad;
            past_subgrads = past_subgrads + its(k-l).subgrad;
        end
        
        z = xk - ak * past_gradients;

        % define objective for prox operator
        if(min(d) == 1) % x = vector
            Obj = past_subgrads' * (x_sdp - z) + ...
                  1/(2*ak) * norm(x_sdp - z, 2)^2;
        else            % x = matrix
            Obj = past_subgrads(:)' * (x_sdp(:) - z(:)) + ...
                  1/(2*ak) * norm(x_sdp - z, 'fro')^2;
        end
        
        info = optimize([options.default_constraint, ...
                                constraints{i}], Obj, options.opt_settings);
        assert(not(info.problem), sprintf('Problem optimizing: %s', info.info));

        % store x(k+1) and its gradient and objective function
        x = value(x_sdp);
        its(k+1).x = x;
        its(k+1).f = f(x);
        its(k+1).subgrad = -past_subgrads + (1/ak)*(z-x);
        its(k+1).grad = gradient_f(x);
        its(k+1).time = toc;
        
        if verbose
            p.ping(); 
        end
    
    end
    
    xstar = x;
end