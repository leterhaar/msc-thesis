function res = svd_rank(A, max_factor)
% returns the number of singular values before a gap that is factor big
if nargin < 2
    max_factor = 1e2;
end
singvals = svd(A);
factors = singvals(1:end-1) ./ singvals(2:end);
res = find(factors > max_factor, 1);
if isempty(res)
    res = length(singvals);
end
% display(sprintf('\nFactor says: %i, difference says: %i', res, find(svd(A) < 1/max_factor,1)-1));