function res = zero_for_nan(A)
% replaces every nan in A with 0
    [m,n] = size(A);
    A_tmp = A(:);
    nans = isnan(A_tmp);
    numnans = sum(nans);
    if numnans > 1
        A_tmp(nans) = 0;
        res = reshape(A_tmp, m, n);
    else
        res = A;
    end
end