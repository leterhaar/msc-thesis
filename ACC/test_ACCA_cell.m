
if not(exist('svm', 'var'))
    test_init;
end

ops = sdpsettings('verbose', 0, 'solver', 'gurobi');
diam = 1;

G = ones(m) - eye(m);
x0 = -10*rand(svm.d,1);

% run the AACC algorithm
time_all = tic;
descriptions = {'ACCA cell x', 'ACCA normal x'};

[xstar_ACCA1, agents_ACCA1] = ACCA_fcn_cell({svm.B}, svm.deltas, @(x) svm.f(x{1}), ...
                           @(x, delta) svm.cons_fcn(x{1}, delta), ...
                           'opt_settings', ops,...
                           'n_agents', m,...
                           'connectivity', G,...
                           'max_its', 30,...
                           'diameter', diam,...
                           'debug', 1,...
                           'verbose', 1,...
                           'x0', {x0},...
                           'stepsize', @(k) 1e4, ...
                           'residuals', @(x,delta) svm.residual_delta(x{1}, delta));
toc(time_all)

% run the AACC algorithm
[xstar_ACCA2, agents_ACCA2] = ACCA_fcn(svm.B, svm.deltas, svm.f, ...
                           svm.cons_fcn, ...
                           'opt_settings', ops,...
                           'n_agents', m,...
                           'connectivity', G,...
                           'max_its', 30,...
                           'diameter', diam,...
                           'debug', 1,...
                           'verbose', 1,...
                           'x0', x0,...
                           'stepsize', @(k) 1e4,...
                           'residuals', svm.residual_delta);
% 
% % %%run the ACC algorithm
% [xstar_ACC, agents_ACC] = ACC(svm.B, svm.delta, svm.deltas, svm.f, ...
%                            svm.cons_delta, ...
%                            'opt_settings', ops,...
%                            'n_agents', m,...
%                            'connectivity', G,...
%                            'diameter', diam,...
%                            'debug', 1,...
%                            'verbose', 1,...
%                            'x0', x0,...
%                            'residuals', svm.residual_delta);
                       
%% calculate convergence and feasibility
K = length(agents_ACCA1(1).iterations);
convergence_ACCA1 = nan(K,m);
feasibility_ACCA1 = nan(K,m);
time_per_iteration_ACCA1 = nan(K,m);
optimal_objective = svm.f(svm.Bstar);
optimizations_run_ACCA1 = zeros(K,1);
no_cons_used_ACCA1 = nan(K,m);
p = progress('Checking constraints', m);
timing_ACCA1 = zeros(K,1);

for i = 1:m
    for k = 1:K
        % calculate difference with centralized objective
        convergence_ACCA1(k,i) = abs(agents_ACCA1(i).iterations(k).J ...
                                                - optimal_objective);
        
        % calculate feasibility percentage
        assign(svm.B, agents_ACCA1(i).iterations(k).x{1})
        feasibility_ACCA1(k,i) = sum(check(svm.cons) < -1e-6) / N_svm * 100;
        
        % store times
        time_per_iteration_ACCA1(k,i) = agents_ACCA1(i).iterations(k).time;
        timing_ACCA1(k) = timing_ACCA1(k) + agents_ACCA1(i).iterations(k).time;

        
        % store total number of iterations run
        if k > 1
            optimizations_run_ACCA1(k) = optimizations_run_ACCA1(k) + ...
                                    agents_ACCA1(i).iterations(k).info.optimized;
            
            no_cons_used_ACCA1(k,i) = agents_ACCA1(i).iterations(k).info.num_cons;

        end
    end
    p.ping();
end

%% calculate convergence and feasibility
K = length(agents_ACCA2(1).iterations);

convergence_ACCA2 = nan(K,m);
feasibility_ACCA2 = nan(K,m);
time_per_iteration_ACCA2 = nan(K,m);
optimal_objective = svm.f(svm.Bstar);
optimizations_run_ACCA2 = zeros(K,1);
no_cons_used_ACCA2 = nan(K,m);
p = progress('Checking constraints', m);
timing_ACCA2 = zeros(K,1);
for i = 1:m
    for k = 1:K
        % calculate difference with centralized objective
        convergence_ACCA2(k,i) = abs(agents_ACCA2(i).iterations(k).J ...
                                                - optimal_objective);
        
        % calculate feasibility percentage
        assign(svm.B, agents_ACCA2(i).iterations(k).x)
        feasibility_ACCA2(k,i) = sum(check(svm.cons) < -1e-6) / N_svm * 100;
        
        % store times
        time_per_iteration_ACCA2(k,i) = agents_ACCA2(i).iterations(k).time;
        timing_ACCA2(k) = timing_ACCA2(k) + agents_ACCA2(i).iterations(k).time;

        % store total number of iterations run
        if k > 1
            optimizations_run_ACCA2(k) = optimizations_run_ACCA2(k) + ...
                                agents_ACCA2(i).iterations(k).info.optimized;
            no_cons_used_ACCA2(k,i) = agents_ACCA2(i).iterations(k).info.num_cons;
        end
    end
    p.ping();
end

%% plot
fig = initfig('ACC iterations', 4);
ax = subplot(211, 'YScale', 'log');
grid on
hold on
hs1 = plot(convergence_ACCA1, 'color', green);
hs2 = plot(convergence_ACCA2, ':', 'color', blue);
legend([hs1(1) hs2(1)], descriptions)
ylabel('|f(x_k^i) - f(x^*) |')
title(sprintf('ACC convergence for SVM with d=%i, m=%i ACCA', d, N_svm));

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
grid on
hold on
plot(feasibility_ACCA1, 'color', green);
plot(feasibility_ACCA2, ':', 'color', blue);

ylabel('% violated');
xlabel('iterations');

initfig('ACC timing', 2);
ax = subplot(211);
grid on
hold on
hs1 = plot(time_per_iteration_ACCA1, 'o', 'color', green);
hs2 = plot(time_per_iteration_ACCA2, '+', 'color', blue);
plot(timing_ACCA1, 'o-', 'color', green);
plot(timing_ACCA2, ':+', 'color', blue);
legend([hs1(1) hs2(1)], descriptions)
ylabel('Time per iteration');

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
grid on
hold on
plot(optimizations_run_ACCA1, 'color', green);
plot(optimizations_run_ACCA2, ':', 'color', blue);
ylabel('Optimizations run');
xlabel('Iteration');

initfig('No of constraints used', 3);
hs1 = plot(no_cons_used_ACCA1, 'color', green);
hs2 = plot(no_cons_used_ACCA2, ':', 'color', blue);
legend([hs1(1) hs2(1)], descriptions)

uistack(fig);
figure(fig);
%% make assertions
verify(all_close(xstar_ACCA1, svm.Bstar, 1e-4), 'ACCA1 not optimal');
verify(all_close(xstar_ACCA2, svm.Bstar, 1e-4), 'ACCA2 not optimal');
verify(sum(feasibility_ACCA1(end,:)) == 0, 'ACCA1 not all feasible');
verify(sum(feasibility_ACCA2(end,:)) == 0, 'ACCA2 not all feasible');