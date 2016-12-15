
ops = sdpsettings('solver', 'mosek');
diam = 2;
G = random_graph(m, diam, 'rand');

% run the AACC algorithm
[xstar_ACCA, agents_ACCA] = ACCA(svm.B, svm.delta, svm.deltas, svm.f, ...
                           svm.cons_delta, ...
                           'opt_settings', ops,...
                           'n_agents', m,...
                           'connectivity', G,...
                           'max_its', 20,...
                           'diameter', diam,...
                           'debug', 1,...
                           'verbose', 1);

% run the ACC algorithm
[xstar_ACC, agents_ACC] = ACC(svm.B, svm.delta, svm.deltas, svm.f, ...
                           svm.cons_delta, ...
                           'opt_settings', ops,...
                           'n_agents', m,...
                           'connectivity', G,...
                           'diameter', diam,...
                           'debug', 1,...
                           'verbose', 1);
                       
%% calculate convergence and feasibility
K = length(agents_ACCA(1).iterations);
convergence_ACCA = nan(K,m);
feasibility_ACCA = nan(K,m);
time_per_iteration_ACCA = nan(K,m);
optimal_objective = svm.f(svm.Bstar);
optimizations_run_ACCA = zeros(K,1);
no_cons_used_ACCA = nan(K,m);
p = progress('Checking constraints', m);
for i = 1:m
    for k = 1:K
        % calculate difference with centralized objective
        convergence_ACCA(k,i) = abs(agents_ACCA(i).iterations(k).J ...
                                                - optimal_objective);
        
        % calculate feasibility percentage
        assign(svm.B, agents_ACCA(i).iterations(k).x)
        feasibility_ACCA(k,i) = sum(check(svm.cons) < -1e-6) / N_svm * 100;
        
        % store times
        time_per_iteration_ACCA(k,i) = agents_ACCA(i).iterations(k).time;
        
        % store total number of iterations run
        if k > 1
            optimizations_run_ACCA(k) = optimizations_run_ACCA(k) + ...
                                    agents_ACCA(i).iterations(k).info.optimized;
            
            no_cons_used_ACCA(k,i) = agents_ACCA(i).iterations(k).info.num_cons;

        end
    end
    p.ping();
end

%% calculate convergence and feasibility
K = length(agents_ACC(1).iterations);

convergence_ACC = nan(K,m);
feasibility_ACC = nan(K,m);
time_per_iteration_ACC = nan(K,m);
optimal_objective = svm.f(svm.Bstar);
optimizations_run_ACC = zeros(K,1);
no_cons_used_ACC = nan(K,m);
p = progress('Checking constraints', m);
for i = 1:m
    for k = 1:K
        % calculate difference with centralized objective
        convergence_ACC(k,i) = abs(agents_ACC(i).iterations(k).J ...
                                                - optimal_objective);
        
        % calculate feasibility percentage
        assign(svm.B, agents_ACC(i).iterations(k).x)
        feasibility_ACC(k,i) = sum(check(svm.cons) < -1e-6) / N_svm * 100;
        
        % store times
        time_per_iteration_ACC(k,i) = agents_ACC(i).iterations(k).time;
        
        % store total number of iterations run
        if k > 1
            optimizations_run_ACC(k) = optimizations_run_ACC(k) + ...
                                agents_ACC(i).iterations(k).info.optimized;
            no_cons_used_ACC(k,i) = agents_ACC(i).iterations(k).info.num_cons;
        end
    end
    p.ping();
end

%% plot
initfig('ACC iterations', 1);
ax = subplot(211, 'YScale', 'log');
grid on
hold on
plot(convergence_ACCA, 'color', green);
plot(convergence_ACC, ':', 'color', blue);
ylabel('|f(x_k^i) - f(x^*) |')
title(sprintf('ACC convergence for SVM with d=%i, m=%i ACCA', d, N_svm));

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
grid on
hold on
plot(feasibility_ACCA, 'color', green);
plot(feasibility_ACC, ':', 'color', blue);
ylabel('% violated');
xlabel('iterations');

initfig('ACC timing', 2);
ax = subplot(211);
grid on
hold on
plot(time_per_iteration_ACCA, 'color', green);
plot(time_per_iteration_ACC, ':', 'color', blue);
ylabel('Time per iteration');

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
grid on
hold on
plot(optimizations_run_ACCA, 'color', green);
plot(optimizations_run_ACC, ':', 'color', blue);
ylabel('Optimizations run');
xlabel('Iteration');

initfig('No of constraints used', 3);
plot(no_cons_used_ACCA, 'color', green);
plot(no_cons_used_ACC, ':', 'color', blue);
%% make assertions
assert(all_close(xstar_ACC, svm.Bstar, 1e-4), 'Not close');
assert(all_close(xstar_ACCA, svm.Bstar, 1e-4), 'Not close');
assert(sum(feasibility_ACC(k,:)) == 0, 'Not all feasible');
assert(sum(feasibility_ACCA(k,:)) == 0, 'Not all feasible');