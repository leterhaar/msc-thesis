function res = normV(vec, p)
% returns the p-normalized version of the vector
[m,n] = size(vec);
assert(n == 1 || m == 1, 'normV only accepts a vector');
if nargin < 2
    p = 1;
end
res = vec ./ norm(vec, p);