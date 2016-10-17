function res = all_close(A, B, tol)
% compares A and B with tolerance
if nargin < 3
    tol = 1e-6;
end
res = all(abs(A(:) - B(:)) < tol);