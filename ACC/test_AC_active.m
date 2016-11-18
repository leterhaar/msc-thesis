tol = 1e-6;
for i = 1:3
    [residuals, Ws] = AC_g(random_x, ac, wind.slice(i), 1);
    params = AC_active(random_x, ac, wind.slice(i), 1);
    all_params = 1:length(residuals);
    active_params = all_params(residuals > -tol & residuals < tol);
    assert(all(active_params == params), 'no match');
    for j = 1:length(residuals)
        param = AC_active(random_x, ac, wind.slice(i), 1, j);
        if residuals(j) > -tol && residuals(j) < tol
            assert(param == j, 'should be passed through');
        else
            assert(isempty(param), 'is not active');
        end
    end
end