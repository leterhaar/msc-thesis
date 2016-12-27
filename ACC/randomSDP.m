%% Test the ACC algorithm for a random SDP
function randomSDP
    yalmip('clear') 

    % define settings
    d = 4;          % dimension of decision variable
    m = 3;          % dimension of constraint set per 1 uncertainty
    N = 100;         % number of scenarios      

    % create SDP vars for every scenario
    x = {sdpvar(d), sdpvar(d), sdpvar(d)};
    delta = sdpvar(m, 1);
    
    % create scenarios
    deltas = rand(N,m);

    % create constraints
    C_all = [x{1} >= 0, x{2} >= 0, x{3} >= 0];
    As = cell(m,1);
    for i = 1:m
        As{i} = randsymmatrix(d);
        for j = 1:N
            C_all = [C_all, trace(As{i}*x{1}) <= deltas(j,i)];
            C_all = [C_all, trace(As{i}*x{2}) <= 2*deltas(j,i)];
            C_all = [C_all, trace(As{i}*x{3}) <= 3*deltas(j,i)];

        end
    end
    cons_fcn_handle = @(x, delta) cons_fcn(x, delta, As);

    
    % check that they are the same
    C_all2 = [x{1} >= 0, x{2} >= 0, x{3} >= 0];
    for i = 1:N
        C_all2 = [C_all2, cons_fcn_handle(x, deltas(i, :))];
    end
    
    assign_cell(x, {randsymmatrix(d), randsymmatrix(d), randsymmatrix(d)});
    res1 = check(C_all);
    res2 = check(C_all2);
    verify(all_close(sort(res1), sort(res2)), 'Not the same');
    
    % make objective
    A0 = randsymmatrix(d);
    obj_fcn = @(x) trace(A0*x{1})+trace(A0*x{2})+trace(A0*x{3});
    Obj = obj_fcn(x);
    
    % solve centralized
    tic
    ops = sdpsettings('verbose', 0, 'solver', 'mosek');
    status = optimize(C_all, Obj, ops);
    assert(not(status.problem), status.info);
    xstar_centralized = values_cell(x);
    fprintf('Centralized solved in %g seconds\n', toc);
    
    % run ACC algorithm
    [xstar, agents] = ACC_fcn_cell(x, deltas, obj_fcn, cons_fcn_handle, ...
                              'verbose', 1,...
                              'debug', 1,...
                              'opt_settings', ops,...
                              'default_constraint', [x{1} >= 0, x{2} >= 0, x{3} >= 0]);
    verify(all_close(xstar_centralized, xstar), 'Not the same');
    check_and_plot(agents, obj_fcn(xstar_centralized), @(x_values) residuals_fcn(x, x_values, C_all));
end

function C = cons_fcn(x, delta, As)
% creates the constraints
    C = [];
    for i = 1:length(delta)
        C = [C, trace(As{i}*x{1}) <= delta(i)];
        C = [C, trace(As{i}*x{2}) <= 2*delta(i)];
        C = [C, trace(As{i}*x{3}) <= 3*delta(i)];
    end
end

function A = randsymmatrix(d)
% create a random symmetric matrix
    A = random_symmetric_matrix(d);
    return
    a = rand(d,1);
    A = a * a';
end

function res = residuals_fcn(x_sdp, x_values, C_all)
% return residuals
    assign_cell(x_sdp, x_values);
    res = check(C_all);
end
