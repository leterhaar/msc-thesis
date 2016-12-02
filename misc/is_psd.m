function res = is_psd(A, tol)
% Returns true if all eigenvalues of A are positive with tolerance
if nargin < 2
    tol = 1e-3;
end
res = all(eig(A) > -tol);