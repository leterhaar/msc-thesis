function A = Ac(C, tol)
% returns the active constraints

if nargin < 2
    tol = 1e-3;
end

% check for infeasibility
assert(all(check(C) > -tol), 'C has infeasible constraints!');

% return tight constraints
A = find((check(C) < tol));