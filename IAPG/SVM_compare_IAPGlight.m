%% Initialize
addpath('../misc');
addpath('../formulation_SVM');
yalmip('clear');
clear
d = 200;                 % dimension of x
m = 1000;                 % number of constraints
max_its = 250;         % maximum number of iterations
alpha = @(k) 1/(k+1);  % step size function

opt_settings = sdpsettings('verbose', 0, 'solver', 'gurobi');

svm = create_SVM(d,m);

% repeat functions and constraints
fs = cell(m,1);
[fs{:}] = deal(@(B) svm.f(B) ./ m);
grad_fs = cell(m,1);
[grad_fs{:}] = deal(@(B) svm.grad_f(B) ./ m);
constraints = cell(m,1);
for i = 1:m
    constraints{i} = svm.cons(i);
end

% x0 = -10*rand(d,1);
% x0 = [];
optimal_objective = svm.f(svm.Bstar);
%% call problem solver
[Bstar_IPG, its_IAPG_light_register] = IAPG_light_register(svm.B, fs, grad_fs, constraints, ...
                            'x0', x0, ...
                            'verbose', 1, ...
                            'max_its', max_its, ...
                            'alpha', alpha,...
                            'opt_settings', opt_settings);

[Bstar_IAPG, its_IAPG_light] = IAPG_light(svm.B, fs, grad_fs, constraints, ...
                            'x0', x0, ...
                            'verbose', 1, ...
                            'max_its', max_its, ...
                            'alpha', alpha, ...
                            'opt_settings', opt_settings);


[Bstar_IAPG2, its_IAPG] = IAPG(svm.B, fs, grad_fs, constraints, ...
                            'x0', x0, ...
                            'verbose', 1, ...
                            'max_its', max_its, ...
                            'alpha', alpha, ...
                            'opt_settings', opt_settings);

%% check feasibility a posteriori
tol = 1e-6;
p = progress('Checking constraints', max_its);
residuals = nan(max_its,1);
for k = 1:max_its
    assign(svm.B, its_IAPG_light_register(k).x);
    residuals = check(svm.cons);
    its_IAPG_light_register(k).feas = sum(residuals < -1e-6) / m * 100;
    
    assign(svm.B, its_IAPG(k).x);
    residuals = check(svm.cons);
    its_IAPG(k).feas = sum(residuals < -1e-6) / m * 100;
    
    assign(svm.B, its_IAPG_light(k).x);
    residuals = check(svm.cons);
    its_IAPG_light(k).feas = sum(residuals < -1e-6) / m * 100;
    
    p.ping();
end
%% calculate differences
differences_IAPG_register = arrayfun(@(i) abs(svm.f(its_IAPG_light_register(i).x) ...
                                          - optimal_objective), 1:max_its);
differences_IAPG = arrayfun(@(i) abs(svm.f(its_IAPG(i).x) ...
                                          - optimal_objective), 1:max_its);
differences_IAPG_light = arrayfun(@(i) abs(svm.f(its_IAPG_light(i).x) ...
                                          - optimal_objective), 1:max_its);
%%
%  plot iterations
initfig('SVM comparison', 2);
hold off
ax = subplot(211);
semilogy(differences_IAPG_register);
hold on
semilogy(differences_IAPG_light, 'linestyle', '-.');
semilogy(differences_IAPG);
grid on;
title(sprintf('SVM Objectives - step-size %s, d=%i, m=%i, infeasible x0, cyclic', ...
                                                    funcname(alpha), d, m));
ylabel('|f(x_k) - f(x^*)|');
legend('IAPG-l-r','IAPG-l','IAPG');

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
hold on
plot([its_IAPG_light_register.feas]');
plot([its_IAPG_light.feas]', 'linestyle', '-.');
plot([its_IAPG.feas]');
grid on;
ylabel('% violated');
title('Feasibility');
%% 
klaarrr