% time the ACC(A) problem for a varying number of constraints and problem
% dimension

yalmip('clear');

dimensions = [20 50 100];
constraint_dimensions = [200 500 1000 2000];
ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
ACCA_times = nan(length(dimensions), length(constraint_dimensions));
ACC_times = nan(length(dimensions), length(constraint_dimensions));
i = 1;
% loop over various problem dimensions
for d = dimensions
    j = 1;
    % loop over various constraint sets
    for N = constraint_dimensions
        
        svm = create_SVM(d,N);
        m = ceil(N/20); % 20 scenarios per agent initially
        G = ones(m) - eye(m);  % fully connected graph
        diam = 1;
        
        % run the ACCA
        tic
        [xstar_ACCA, agents_ACCA] = ACCA_fcn(svm.B, svm.deltas, svm.f, ...
                           svm.cons_fcn, ...
                           'opt_settings', ops,...
                           'n_agents', m,...
                           'connectivity', G,...
                           'max_its', 50,...
                           'diameter', diam,...
                           'debug', 1,...
                           'verbose', 0,...
                           'residuals', svm.residual_delta);
        ACCA_times(i,j) = 1000*toc;
        
        % check optimiality and feasibility
        verify(all_close(xstar_ACCA, svm.Bstar, 1e-3), sprintf('ACCA d%i N%i not optimal', d, N), 0);
        assign(svm.B, xstar_ACCA);
        residuals = check(svm.cons);
        verify(sum(residuals < -1e-6) == 0, sprintf('ACCA d%i N%i not feasible', d, N), 0);

        %%run the ACC algorithm
        tic
        [xstar_ACC, agents_ACC] = ACC(svm.B, svm.delta, svm.deltas, svm.f, ...
                           svm.cons_delta, ...
                           'opt_settings', ops,...
                           'n_agents', m,...
                           'connectivity', G,...
                           'diameter', diam,...
                           'debug', 1,...
                           'verbose', 0,...
                           'residuals', svm.residual_delta);
        ACC_times(i,j) = 1000*toc;
        
        % check optimiality and feasibility
        verify(all_close(xstar_ACC, svm.Bstar, 1e-3), sprintf('ACC d%i N%i not optimal', d, N), 0);
        assign(svm.B, xstar_ACC);
        residuals = check(svm.cons);
        verify(sum(residuals < -1e-6) == 0, sprintf('ACC d%i N%i not feasible', d, N), 0);

        j = j + 1;
    end
    i = i + 1;
end

%% plot times
h = initfig('Opt times vs problem dimensions ACCA/ACC', 1);
subplot(211)
bar3(ACC_times);
title('ACC');
ylabel('Problem dimension');
xlabel('Constraint set size');
set(gca, 'xticklabels', constraint_dimensions);
set(gca, 'yticklabels', dimensions);
zlabel('Time to solve (ms)');

subplot(212)
bar3(ACCA_times);
title('ACCA');
ylabel('Problem dimension');
xlabel('Constraint set size');
set(gca, 'xticklabels', constraint_dimensions);
set(gca, 'yticklabels', dimensions);
zlabel('Time to solve (ms)');
uistack(h)
