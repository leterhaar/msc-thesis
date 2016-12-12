addpath('../formulation_SVM');

d = 50;
m = 100;

svm = create_SVM(d,m);

% repeat functions and constraints
fs = cell(m,1);
[fs{:}] = deal(@(x) svm.f(x) ./ m);
grad_fs = cell(m,1);
[grad_fs{:}] = deal(@(x) svm.grad_f(x) ./ m);
constraints = cell(m,1);
for i = 1:m
    constraints{i} = svm.cons(i);
end

% start with a bad x0
x0 = -10*rand(d,1);


% algorithm settings
opt_settings = sdpsettings('verbose', 0, 'solver', 'gurobi');
max_its = 1000;
alpha = @(k) 1/(k+10);
[Bstar, its] = IAPG(svm.B, fs, grad_fs, constraints, ...
                    'max_its', max_its, ...
                    'alpha', alpha, ...
                    'opt_settings', opt_settings, ...
                    'verbose', 1,...
                    'x0', x0);
%% plot iterations
optimal_objective = svm.f(svm.Bstar);
differences = abs((m.*[its.f])-optimal_objective);

initfig('SVM IAPG Objective', 2);
ax = subplot(211);
semilogy(differences);
hold on;
grid on;
title('SVM IAPG Objective');
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

initfig('Maximum delay');
plot([its.b]');

                                        