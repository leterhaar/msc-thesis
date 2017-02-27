%% [xstar, iterations] = PG(x_sdp, F, grad_F, constraints, ...)
% runs the incremtal aggregated proximal gradient algorithm and returns 
% solution and iterations
% 
% solves min F(x) s.t. x \in X for i = 1,...,m
%
% PARAMETERS
% * x               : an sdpvar
% * F               : a single function handle F(x)
% * grad_F          : a cell with m function handle to the gradient of the smooth part
%                     of the objective function gradient_f(x)
% * constraints     : all the constraints
% * options         : set of options, optional, as follows
%   - x0            : value for x0. if empty, a feasible x0 will be found
%   - verbose       : 0 (hide) or 1 (show progress, default)
%   - max_its       : maximum number of iterations (default 500)
%   - opt_settings  : optimization settings (default verbose = 0)
%   - alpha         : function handle for alpha. default alpha = 1/(k+1)

function [xstar, its] = ProxGrad(x_sdp, F, grad_F, constraints, varargin)

    % check types of input
    assert(isa(x_sdp, 'sdpvar'), 'x should be sdpvar');
    assert(isa(F, 'function_handle'), 'f should be function handle');
    assert(isa(grad_F, 'function_handle'), ...
                                   'gradient_f should be function handle');
    assert(isa(constraints, 'constraint') || isa(constraints, 'lmi'), ...
                                      'constraints should be constraints');
    
    % define default options
    options = struct('verbose', 1, 'max_its', 500, 'opt_settings', ...
                     sdpsettings('verbose', 0), 'x0', [], ...
                     'default_constraint', [], 'batch', 1, ...
                     'alpha', @(k) 1/(k+1));
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
    
    %% initialize algorithm
    
    % find feasible x0
    if isempty(options.x0)
        info = optimize(constraints, [], options.opt_settings);
        assert(not(info.problem), 'Problem optimizing: %s', info.info);
        x0 = value(x_sdp);
    else
        assert(all(size(options.x0) == d), 'x0 has the wrong size');
        x0 = options.x0;
    end
    
    its = struct(   'x', [], ...        value for x
                    'z', [], ...        value for z
                    'f', [], ...        objective function value
                    'grad', [], ...     gradient at x
                    'time', []...       time per iterations    
                    );
                
    its(1).x = x0;
    its(1).f = F(x0);
    its(1).grad = grad_F(x0);
    its(1).time = nan;
    
    % cache optimizer
    alpha_sdp = sdpvar(1);
    z_sdp = sdpvar(d(1), d(2), 'full');
    
    if isempty(options.opt_settings.solver)
        warning('Calling optimizer with no solver specified!');
    end
    
    proximal = optimizer(constraints, ...
                         1/(2*alpha_sdp) * norm(x_sdp(:) - z_sdp(:))^2,...
                         options.opt_settings, {alpha_sdp, z_sdp}, x_sdp);
    
    %% start iterating
    if verbose
        p = progress('Iterating PG', options.max_its);
    end
    
    for k = 1:options.max_its
        
        tic
        % gradient step
        its(k).z = its(k).x - options.alpha(k) * its(k).grad;
        
        % projection step
        [x, problem, msg] = proximal(options.alpha(k), its(k).z);
        assert(not(problem), sprintf('Problem optimizing: %s', msg{:}));
        
        its(k+1).x = x;
        its(k+1).f = F(x);
        its(k+1).grad = grad_F(x);
        its(k+1).time = toc;
        
        if verbose
            p.ping();
        end
    end
    
    xstar = x;
end
