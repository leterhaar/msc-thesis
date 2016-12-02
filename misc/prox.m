function next_x = prox(f, alpha, z, x, cons)
% x = prox(f, alpha, z, x, cons)
% evaluates the proximal operator for point z over f with 
% stepsize alpha w.r.t. cons (optional)

    if nargin < 5
        cons = [];
    end
    
    assert(isa(f, 'sdpvar'), 'F should be a SDP var, not a %s', class(f));    
    assert(isa(x, 'sdpvar'), 'x should be a SDP var, not a %s', class(x));
    assert(all(size(x) == size(z)), 'Should be of same dimensions');
    
    if min(size(x)) == 1
        % vector valued problem
        f_prox = f + 1/(2*alpha) * norm(x - z)^2;
    else
        % matric valued problem
        f_prox = f + 1/(2*alpha) * norm(x - z, 'fro')^2;
    end
    
    diagnostics = optimize(cons, f_prox, sdpsettings('verbose', 0));
    assert(not(diagnostics.problem), 'Problem optimizing: %s', diagnostics.info);

    next_x = value(x);
end