function res = is_scalar(A)
% bool = is_scalar(A)
% return true if size of A is [1 1]
    res = all(size(A) == [1 1]);
end