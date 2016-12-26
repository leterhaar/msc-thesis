% time the ACC(A) problem for a varying number of constraints and problem
% dimension

yalmip('clear');
if not(exist('create_SVM', 'file'))
    addpath('../formulation_SVM');
    addpath('../misc');
end


dimensions = [100];
constraint_dimensions = [1000 2000 5000 10000];
ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
tol = 1e-6;
ACCA_times = nan(length(dimensions), length(constraint_dimensions));
ACC_times = nan(length(dimensions), length(constraint_dimensions));
i = 1;
% loop over various problem dimensions
for d = dimensions
    j = 1;
    % loop over various constraint sets
    for N = constraint_dimensions
        
        svm = create_SVM(d,N);
        m = ceil(N/100);        % 100 scenarios per agent initially
        G = ones(m) - eye(m);   % fully connected graph
        diam = 1;
        
        % run the ACCA
        ACCA_start = tic;
        [xstar_ACCA, agents_ACCA] = ACCA_fcn2(svm.B, svm.deltas, svm.f, ...
                           svm.cons_fcn, ...
                           'opt_settings', ops,...
                           'n_agents', m,...
                           'connectivity', G,...
                           'diameter', diam,...
                           'debug', 1,...
                           'verbose', 0,...
                           'residuals', svm.residual_delta);
        ACCA_times(i,j) = 1000*toc(ACCA_start);
        
        % check optimiality and feasibility
        verify(all_close(xstar_ACCA, svm.Bstar, 1e-5), 'ACCA d%i N%i not optimal: %g', d, N, norm(svm.Bstar-xstar_ACCA, 'inf'));
        verify(all(svm.residuals(xstar_ACCA) > -tol), 'ACCA d%i N%i not feasible', d, N);

        %%run the ACC algorithm
        ACC_start = tic;
        [xstar_ACC, agents_ACC] = ACC(svm.B, svm.delta, svm.deltas, svm.f, ...
                           svm.cons_delta, ...
                           'opt_settings', ops,...
                           'n_agents', m,...
                           'connectivity', G,...
                           'diameter', diam,...
                           'debug', 1,...
                           'verbose', 0,...
                           'residuals', svm.residual_delta);
        ACC_times(i,j) = 1000*toc(ACC_start);
        
        % check optimiality and feasibility
        verify(all_close(xstar_ACC, svm.Bstar, 1e-5), 'ACC d%i N%i not optimal: %g', d, N, norm(svm.Bstar-xstar_ACCA, 'inf'));
        verify(all(svm.residuals(xstar_ACC) > -tol), 'ACC d%i N%i not feasible', d, N);
        j = j + 1;
    end
    i = i + 1;
end

%% plot times
h = initfig('Opt times vs problem dimensions ACCA/ACC', 1);
d_plot = 100;
N_plot = 1000;
factor = 6e4;       % to make minutes
% subplot(211)
% box on; hold on; grid on;
% plot(dimensions, ACC_times(:, constraint_dimensions == N_plot)./factor);
% plot(dimensions, ACCA_times(:, constraint_dimensions == N_plot)./factor);
% title(sprintf('Solving time vs problem dimension ACC/ACCA (N=%i)', N_plot));
% ylabel('Solving time [min]');
% xlabel('Problem dimension d');
% legend('ACC', 'ACCA', 'location', 'nw');
% 
% subplot(212)
% hold off; 
plot(constraint_dimensions, ACC_times(dimensions == d_plot, :)./factor, 'x-');
hold on; grid on; box on; 
plot(constraint_dimensions, ACCA_times(dimensions == d_plot, :)./factor, 'o-');
title(sprintf('Solving time vs constraint dimension ACC/ACCA (d=%i)', d_plot));
ylabel('Solving time [min]');
xlabel('Constraint set size N');
legend('ACC', 'ACCA', 'location', 'nw');

uistack(h)
