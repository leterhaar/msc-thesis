d = 100;                    % dimension of problem
m = 200;                    % number of constraints
alpha = @(k) 1/(k+1);      % stepsize function
max_its = 500;             % maximum nr of iterations
B = 10;                      % batch size

opt_settings = sdpsettings('verbose', 0, 'solver', 'gurobi');

svm = create_SVM(d,m);
x0 = -10*rand(d,1); %zeros(d,1);    % start with infeasible point
                
[Bstar_IPG, its_IPG] = IPG(svm.B, svm.fs, svm.grad_fs, svm.constraints, ...
                          'verbose', 1,...
                          'alpha', alpha, ...
                          'max_its', max_its, ...
                          'opt_settings', opt_settings,...
                          'x0', x0);
                      
[Bstar_IAPGl, its_IAPGl] = IAPG_light(svm.B, svm.fs, svm.grad_fs, svm.constraints, ...
                          'verbose', 1,...
                          'alpha', alpha, ...
                          'max_its', max_its, ...
                          'opt_settings', opt_settings,...
                          'batch_size', B,...
                          'x0', x0);
%% check feasibility a posteriori
tol = 1e-6;
p = progress('Checking constraints', max_its);
residuals = nan(max_its,1);
for k = 1:max_its
    assign(svm.B, its_IAPGl(k).x);
    residuals = check(svm.cons);
    its_IAPGl(k).feas = sum(residuals < -1e-6) / m * 100;
    
    assign(svm.B, its_IPG(k).x);
    residuals = check(svm.cons);
    its_IPG(k).feas = sum(residuals < -1e-6) / m * 100;
    
    p.ping();
end
%% calculate differences
differences_IAPG2 = arrayfun(@(i) abs(svm.f(its_IAPGl(i).x) ...
                                          - optimal_objective), 1:max_its);
differences_IAPG = arrayfun(@(i) abs(svm.f(its_IPG(i).x) ...
                                          - optimal_objective), 1:max_its);
%% plot iterations
initfig('SVM batch comparison', 8);
hold off
ax = subplot(211);
semilogy(differences_IAPG, 'color', blue);
hold on
semilogy(differences_IAPG2, 'color', orange, 'linestyle', '-.');
grid on;
title(sprintf('SVM Objectives - step-size %s, d=%i, m=%i', ...
    funcname(alpha), d, m));
ylabel('|f(x_k) - f(x^*)|');
legend('IPG', 'IAPG light');

ax2 = subplot(212);
linkaxes([ax ax2], 'x');
hold on
plot([its_IPG.feas]', 'color', blue);
plot([its_IAPGl.feas]', 'color', orange, 'linestyle', '-.');
grid on;
ylabel('% violated');
title('Feasibility');
xlabel('Iterations');
ylim([0 25])
%%
klaarrr