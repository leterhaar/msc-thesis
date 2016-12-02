%% [xstar, iterations] = IPG(x, f, gradient_f, constraints, ...)
% runs the incremtal proximal gradient algorithm and returns iterations
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

function [xstar, its] = IPG(x_sdp, f, gradient_f, constraints, varargin)

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
                     sdpsettings('verbose', 0), 'x0', []);
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
        info = optimize([constraints{:}], [], options.opt_settings);
        assert(not(info.problem), 'Problem optimizing: %s', info.info);
        x0 = value(x_sdp);
    else
        assert(all(size(options.x0) == d), 'x0 has the wrong size');
        x0 = options.x0;
    end

    its = struct(   'x', [], ...        value for x
                    'f', [], ...        objective function value
                    'grad', [], ...     gradient at x
                    'time', []...       time per iterations    
                    );

    its(1).x = x0;
    its(1).grad = gradient_f(x0);
    its(1).f = f(x0);
    its(1).time = nan;

    if verbose
        p = progress('Iterating IPG', options.max_its);
    end

    for k = 1:options.max_its
        tic
        xk = its(k).x;
        ak = 0.005;
        i = randi(min(k,m));

        % gradient step at x(k)
        z = xk - ak * its(k).grad;

        % define objective for prox operator
        if(min(d) == 1) % x = vector
            Obj = 1/(2*ak) * norm(x_sdp - z, 2)^2;
        else            % x = matrix
            Obj = 1/(2*ak) * norm(x_sdp - z, 'fro')^2;
        end
        
        % subgradient at x(k+1) / proximal step
        info = optimize(constraints{i}, Obj, options.opt_settings);
        assert(not(info.problem), sprintf('Problem optimizing: %s', info.info));

        % store x(k+1) and its gradient and objective function
        x = value(x_sdp);
        its(k+1).x = x;
        its(k+1).f = f(x);
        its(k+1).grad = gradient_f(x);
        its(k+1).time = toc;
        
        if verbose
            p.ping(); 
        end

    end
    
    xstar = x;
end