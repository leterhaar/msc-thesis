%% [xstar, iterations] = IPG(x, fs, grad_fs, constraints, ...)
% runs the incremtal proximal gradient algorithm and returns iterations
% 
% solves min_ F(x) = sum f_i(x) s.t. x \in X_i for i = 1,...,m
%
% PARAMETERS
% * x               : an sdpvar
% * fs              : a cell of m function handles to the smooth part of 
%                     the objective function f_i(x)
% * grad_fs         : a cell of m function handles to the gradient of the smooth part
%                     of the objective function gradient_f(x)
% * constraints     : a cell of m constraint sets
% * options         : set of options, optional, as follows
%   - x0            : value for x0. if empty, a feasible x0 will be found
%   - verbose       : 0 (default) or 1 (show progress)
%   - max_its       : maximum number of iterations (default 500)
%   - opt_settings  : optimization settings (default verbose = 0)
%   - default_constraint : default constraint for every optimization step
%   - alpha         : function handle for alpha. default alpha = 1/(k+1)



function [xstar, its] = IPG(x_sdp, fs, grad_fs, constraints, varargin)

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
        info = optimize([options.default_constraint; constraints{:}], [],...
                                                    options.opt_settings);
        assert(not(info.problem), 'Problem optimizing: %s', info.info);
        x0 = value(x_sdp);
    else
        assert(all(size(options.x0) == d), 'x0 has the wrong size');
        x0 = options.x0;
    end

    its = struct(   'x', [], ...        value for x
                    'f', [], ...        objective function value
                    'time', []...       time per iterations    
                    );
    i = randi(m);
    its(1).x = x0;
    its(1).f = fs{i}(x0);
    its(1).time = nan;
    
    %% initialize solvers
    
    % create sdpvar for z, ak
    z_sdp = sdpvar(d(1), d(2), 'full');
    a_sdp = sdpvar(1);
    
    % define objective with vectorized versions (so matrix also works)
    Obj = 1/(2*a_sdp) * norm(x_sdp(:) - z_sdp(:), 2)^2;
    proximals = cell(m,1);
    
    if isempty(options.opt_settings.solver)
        warning('Calling optimizer with no solver specified!');
    end
    
    if verbose
        p = progress('Preparing IPG', m);
    end
    for i = 1:m
        proximals{i} = optimizer([options.default_constraint; ...
                                  constraints{i}], Obj, options.opt_settings, ...
                                  {a_sdp, z_sdp}, ...
                                  x_sdp);
        if verbose
            p.ping();
        end
    end


    %% iterate
    if verbose
        p = progress('Iterating IPG', options.max_its);
    end

    for k = 1:options.max_its
        tic
        xk = its(k).x;
        ak = options.alpha(k);
        i = randi(m);

        % gradient step at x(k)
        zk = xk - ak * grad_fs{i}(xk);

        % solve
        [x, problem, msg] = proximals{i}(ak, zk);
        assert(not(problem), sprintf('Problem optimizing: %s', msg{:}));

        % store x(k+1) and its gradient and objective function
        its(k+1).x = x;
        its(k+1).f = fs{i}(x);
        its(k+1).time = toc;
        
        if verbose
            p.ping(); 
        end

    end
    
    xstar = x;
end
