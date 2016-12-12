%% Initialize
addpath('../misc');
addpath('../formulation_SVM');
yalmip('clear');
clear
d = 50;                 % dimension of x
m = 100;                 % number of constraints
max_its = 500;         % maximum number of iterations
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

x0 = -10*rand(d,1);
% x0 = [];
optimal_objective = svm.f(svm.Bstar);
%% call problem solver
[Bstar_IPG, its_IPG] = IPG(svm.B, fs, grad_fs, constraints, ...
                            'x0', x0, ...
                            'verbose', 1, ...
                            'max_its', max_its, ...
                            'opt_settings', opt_settings);

[Bstar_IAPG, its_IAPG_light] = IAPG_light(svm.B, fs, grad_fs, constraints, ...
                            'x0', x0, ...
                            'verbose', 1, ...
                            'max_its', max_its, ...
                            'alpha', alpha, ...
                            'opt_settings', opt_settings);


[Bstar_IAPG2, its_IAPG] = IAPG2(svm.B, fs, grad_fs, constraints, ...
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
    assign(svm.B, its_IPG(k).x);
    residuals = check(svm.cons);
    its_IPG(k).feas = sum(residuals < -1e-6) / m * 100;
    
    assign(svm.B, its_IAPG(k).x);
    residuals = check(svm.cons);
    its_IAPG(k).feas = sum(residuals < -1e-6) / m * 100;
    
    assign(svm.B, its_IAPG_light(k).x);
    residuals = check(svm.cons);
    its_IAPG_light(k).feas = sum(residuals < -1e-6) / m * 100;
    
    p.ping();
end
%% calculate differences
differences_IPG = arrayfun(@(i) abs(svm.f(its_IPG(i).x) ...
                                          - optimal_objective), 1:max_its);
differences_IAPG = arrayfun(@(i) abs(svm.f(its_IAPG(i).x) ...
                                          - optimal_objective), 1:max_its);
differences_IAPG_light = arrayfun(@(i) abs(svm.f(its_IAPG_light(i).x) ...
                                          - optimal_objective), 1:max_its);

%% calculate norm of h_i(x_k+1)
norm_his_IAPG = arrayfun(@(i) norm(its_IAPG(i).subgrad), 1:max_its);
norm_his_IAPG_light = arrayfun(@(i) norm(its_IAPG_light(i).subgrad), 1:max_its);
running_sum_IAPG = ...
               arrayfun(@(i) sum(norm_his_IAPG(max(1,i-m):i)), 1:max_its);
running_sum_IAPG_light = ...
               arrayfun(@(i) sum(norm_his_IAPG_light(max(1,i-m):i)), 1:max_its);
           
norm_his_IAPG_light(norm_his_IAPG_light < 1e-6) = nan;
norm_his_IAPG(norm_his_IAPG < 1e-6) = nan;

%  plot iterations
initfig('SVM comparison', 7);
hold off
ax = subplot(311);
semilogy(differences_IPG);
hold on
semilogy(differences_IAPG_light, 'linestyle', '-.');
semilogy(differences_IAPG);
grid on;
title(sprintf('SVM Objectives - step-size %s, d=%i, m=%i, infeasible x0, cyclic', ...
                                                    funcname(alpha), d, m));
ylabel('|f(x_k) - f(x^*)|');
legend('IPG','IAPG-light','IAPG');

ax2 = subplot(313);

hold on
plot([its_IPG.feas]');
plot([its_IAPG_light.feas]', 'linestyle', '-.');
plot([its_IAPG.feas]');
grid on;
ylabel('% violated');
title('Feasibility');

ylim([0 100])

ax3 = subplot(312);
linkaxes([ax ax2 ax3], 'x');
hold on

h1 = plot(0,0, ':');
h2 = plot(running_sum_IAPG_light, ':');
h3 = plot(running_sum_IAPG, ':');
plot(norm_his_IAPG_light, 'x', 'color', get(h2, 'color'));
plot(norm_his_IAPG, 'o', 'color', get(h3, 'color'));

h = legend( '$\| \tilde \nabla h_i(x_{k+1}) \|$ - IPG', ...
            '$\| \tilde \nabla h_i(x_{k+1}) \|$ - IAPG light', ...
            '$\| \tilde \nabla h_i(x_{k+1}) \|$ - IAPG', ...
            '$ \sum_{j \neq i} \| \tilde \nabla h_j(x_{k+1}) \| $ - IPG', ...
            '$ \sum_{j \neq i} \| \tilde \nabla h_j(x_{k+1}) \| $ - IAPG light', ...
            '$ \sum_{j \neq i} \| \tilde \nabla h_j(x_{k+1}) \| $ - IAPG');
set(h, 'interpreter', 'latex')
grid on
xlabel('Iterations');
title('Magnitude of $\tilde \nabla h_i(x_{k+1})$', 'interpreter', 'latex');
ylabel('$\| \tilde \nabla h_i(x_{k+1}) \| $', 'interpreter', 'latex');
xlim([0 500]);
%%  plot delay
initfig('Max delay', 6);
plot([its_IAPG.b], 'color', blue);
plot([its_IAPG_light.b], 'color', orange, 'linestyle', '-.');
legend('IAPG - start with IPG', 'IAPG - start with aggregation');

%% 
klaarrr