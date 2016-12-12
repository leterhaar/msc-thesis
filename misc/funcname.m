function res = funcname(handle)
% returns function name without input arguments for anonymous functions
    res = func2str(handle);
    if strcmp('@', res(1))
        haakjes = find(res == ')');
        res = res(haakjes(1)+1:end);
    end
end
