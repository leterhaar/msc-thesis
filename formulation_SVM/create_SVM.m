function prob = create_SVM(d, m)
% creates a feasible SVM problem definition 
%
% min 0.5 || B ||^2
% s.t. y_i (x^t * B) >= 1   for i = 1, ..., m
%
% PARAMETERS
% ==========
% d : dimension of x and B
% m : number of data points
%
% RETURNS
% =======
% A structure with the following fields
% xs        : array of d x m data points
% ys        : array of m labels {-1, 1}
% B         : sdpvar(d,1,'full');
% cons      : lmi with all the constraints
% f         : function handle for the objective function: f(B)
% grad_f    : function handle for the gradient function: grad_f(B)
% residual  : function handle for the residuals of the form residual(B,i)
% Bstar     : optimal value for B  
% constraints : cell with single constraints as entries
% delta     : sdpvar for the uncertainty (d(1) = y, d(2:end) = x)
% cons_delta : the constraints as a function of B and delta

    tic
    % loop so that we always have a feasible problem
    while 1
        xs = randn(d, m);
        separation = 3*rand(d,1);
        xs = xs + [repmat(separation, 1, floor(m/2)) ...
                  -repmat(separation, 1, ceil(m/2))];
        ys = [ones(floor(m/2), 1); -ones(ceil(m/2), 1)];

        shuffled_ids = randperm(m);
        prob.xs = xs(:, shuffled_ids);
        prob.ys = ys(shuffled_ids);
        prob.B = sdpvar(d,1,'full');
        prob.delta = sdpvar(1, d+1);
        prob.d = d;
        prob.m = m;

        % define constraints
        prob.cons = [];
        for i = 1:m
            prob.cons = [prob.cons, prob.ys(i) * prob.xs(:, i)' * prob.B >= 1];
        end
        
        % define constraints as function of delta
        prob.cons_delta = prob.delta(1) * prob.delta(2:end) * prob.B >= 1;
        prob.deltas = [prob.ys prob.xs'];

        prob.f = @(B) 0.5*norm(B)^2;
        prob.grad_f = @(B) B;
        prob.residual = @(B,i) prob.ys(i) * prob.xs(:, i)' * B - 1;
        prob.residuals = @(B) arrayfun(@(i) prob.ys(i) * prob.xs(:, i)' * B - 1, 1:m);
        prob.residual_delta = @(B,delta) delta(1) * delta(2:end) * B - 1;

        
        % solve problem and store optimal value
        diagnostics = optimize(prob.cons, prob.f(prob.B), ...
                                             sdpsettings('verbose', 0));
        if not(diagnostics.problem)
            
            prob.Bstar = value(prob.B);
            % repeat functions and constraints
            prob.fs = cell(m,1);
            [prob.fs{:}] = deal(@(B) prob.f(B)./m);
            prob.grad_fs = cell(m,1);
            [prob.grad_fs{:}] = deal(@(B) prob.grad_f(B)./m);
            prob.constraints = cell(m,1);
            for i = 1:m
                prob.constraints{i} = prob.cons(i);
            end

            break
        end
        
        if toc > 30
            error('Took more than 30 seconds to create a SVM problem, aborting');
        end
    end
end