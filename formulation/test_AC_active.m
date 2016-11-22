tol = 1e-6;

%% check if it works with random parameters

% loop over scenarios
for i = 1:N
    
    % get residuals and test activeness (without j)
    [residuals, Ws] = AC_g(random_x, ac, wind.slice(i), t);
    params = AC_active(random_x, ac, wind.slice(i), t);
    all_params = [1:length(residuals)]';
    active_params = all_params(residuals > -tol & residuals < tol);
    
    % add PSDness of Ws to params check
    if min(eig(Ws)) < tol && min(eig(Ws)) > -tol
        active_params = [active_params; -1];
    end
    
    assert(all(active_params == params), 'no match');
    
    % loop over individually and check again
    for j = 1:length(residuals)
        param = AC_active(random_x, ac, wind.slice(i), t, j);
        if residuals(j) > -tol && residuals(j) < tol
            assert(param == j, 'should be passed through');
        else
            assert(isempty(param), 'is not active');
        end
    end
    
    % check PSD constraints
    if min(eig(Ws)) < tol && min(eig(Ws)) > -tol
        Ws_active = 1;
    else
        Ws_active = 0;
    end
    
    Ws_params = AC_active(random_x, ac, wind.slice(i), t, -1);
    assert(not(isempty(Ws_params)) == Ws_active, 'Ws should match');
end

%% check if it works for optimal solution xstar
all_active_params = [];
for i = 1:N
    
    % get residuals and test activeness (without j)
    [residuals, Ws] = AC_g(xstar, ac, wind.slice(i), t);
    params = AC_active(xstar, ac, wind.slice(i), 1);
    all_params = [1:length(residuals)]';
    active_params = all_params(residuals > -tol & residuals < tol);
    
    % add PSDness of Ws to params check
    Ws_active = 0;
    if min(eig(Ws)) < tol && min(eig(Ws)) > -tol
        Ws_active = 1;
        active_params = [active_params; -1];
    end
    
    % append to all
    all_active_params = [all_active_params; ...
                         [ones(length(active_params),1)*i active_params]];
    
    assert(all(active_params == params), 'no match');
    
    % loop over individually and check again
    for j = 1:length(residuals)
        param = AC_active(xstar, ac, wind.slice(i), t, j);
        if residuals(j) > -tol && residuals(j) < tol
            assert(param == j, 'should be passed through');
        else
            assert(isempty(param), 'is not active');
        end
    end
    
    % check Ws individually
    Ws_params = AC_active(xstar, ac, wind.slice(i), t, -1);
    assert(not(isempty(Ws_params)) == Ws_active, 'Should match');
end

%% check if solution with the active params still is the same

% formulate the constraints again
C2 = AC_cons_det(x, ac, wind, t);

for p = all_active_params'
    i = p(1);
    j = p(2);
    
    C2 = [C2, AC_cons_scen(x, ac, wind.slice(i), t, j)];
end

% optimize and extract solution
diagnostics2 = optimize(C2, Obj, sdpsettings('verbose', 0));
assert(not(diagnostics2.problem), 'Problem optimizing');
xstar2 = values_cell(x);

% check if solutions are the same
assert(all_close(xstar, xstar2), 'Should have same solution');

% check if the solution is feasible and makes sense
problem = AC_check(xstar2, ac, wind, t);
assert(not(problem), sprintf('Problem should be feasible, error %i', problem));


