% test the PG solver with a SVM problem

x0 = -10*rand(d,1);
max_its = 500;

opt_settings = sdpsettings('verbose', 0, 'solver', 'gurobi');

[Bstar, its] = PG(svm.B, svm.f, svm.grad_f, svm.cons, 'max_its', max_its, ...
                   'x0', x0, 'alpha', @(k) 1/(k+10), ...
                   'opt_settings', opt_settings);
% assert(all_close(svm.Bstar, Bstar), 'Does not converge');
%% plot iterations
optimal_objective = svm.f(svm.Bstar);
differences = abs([its.f]-optimal_objective);

initfig('SVM PG Objective', 1);
ax = subplot(211);
semilogy(differences);
hold on;
grid on;
title('SVM PG Objective');
ylabel('|f(x_k) - f(x^*)|');


% check constraints and plot constraint violation
residuals = nan(m, 1);
percentage_violated = nan(1, max_its);
p = progress('Checking constraints', max_its*m);
for k = 1:max_its
    for i = 1:m
        residuals(i) = svm.residual(its(k).x, i);
        p.ping();
    end
    percentage_violated(k) = sum(residuals < -1e-6) / m * 100;
end

ax2 = subplot(212);

linkaxes([ax ax2], 'x');
plot(percentage_violated);
hold on;
grid on;
ylabel('% violated');
title('Feasibility');
xlabel('Iterations');
