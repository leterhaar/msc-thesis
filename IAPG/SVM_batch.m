%% Initialize
addpath('../misc');
yalmip('clear');
clear
d = 25;         % dimension of x
N = 100;        % number of constraints
max_its = 4000; % maximum number of iterations
b = N;          % upper bound on delay 
xs = randn(N, d);
mean1 = 3*rand(1,d);
xs = xs + [repmat(mean1, N/2, 1); -repmat(mean1, N/2, 1)];
ys = [ones(N/2, 1); -ones(N/2, 1)];

%% Create problem

B = sdpvar(d,1, 'full');
Obj = 0.5*norm(B,2)^2;
C_cent = [];
constraints = cell(N,1);
for i = 1:N
    C_cent = [C_cent, ys(i) * xs(i, :) * B >= 1];
    constraints{i} = ys(i) * xs(i, :) * B >= 1;
end
tic
opt_settings = sdpsettings('verbose', 0);
diagnostics = optimize(C_cent, Obj, opt_settings);
assert(not(diagnostics.problem), diagnostics.info);
Bstar_cent = value(B);
toc

% define B0, somewhere in the feasible set
diagnostics = optimize(C_cent, [], opt_settings);
assert(not(diagnostics.problem));
B0 = value(B);

%% call problem solver

% [Bstar_IPG, its_IPG] = IPG(B, @SVM_f, @SVM_gradient_f, constraints, ...
%                             'x0', B0, ...
%                             'verbose', 1, ...
%                             'max_its', max_its);

[Bstar_IAPG, its_IAPG] = IAPG(B, @SVM_f, @SVM_gradient_f, constraints, ...
                            'x0', B0, ...
                            'verbose', 1, ...
                            'max_its', max_its);
% 
%                           
% [Bstar_IAPG2, its_IAPG2] = IAPG_light(B, @SVM_f, @SVM_gradient_f, constraints, ...
%                            'x0', B0, ...
%                            'verbose', 1, ...
%                            'max_its', max_its, ...
%                            'b', b);

%% check feasibility a posteriori
tol = 1e-6;
p = progress('Checking constraints', max_its);
for i = 1:max_its
    
%     % check IPG
%     assign(B, its_IPG(i).x);
%     residuals = check(C_cent);
%     its_IPG(i).feas = sum(residuals < -tol)/N * 100;
    
    % check IAPG
    assign(B, its_IAPG(i).x);
    residuals = check(C_cent);
    its_IAPG(i).feas = sum(residuals < -tol)/N * 100;
%     
%     % check IAPG2
%     assign(B, its_IAPG2(i).x);
%     residuals = check(C_cent);
%     its_IAPG2(i).feas = sum(residuals < -tol)/N * 100;
    
    p.ping();
end
%% plot objectives
initfig('Objectives', 1);
% differences_objective_IPG = abs([its_IPG(2:end).f] - SVM_f(Bstar_cent));
differences_objective_IAPG = abs([its_IAPG(2:end).f] - SVM_f(Bstar_cent));
% differences_objective_IAPG2 = abs([its_IAPG2(2:end).f] - SVM_f(Bstar_cent));
differences_objective_IAPG(1:N) = nan;
% differences_objective_IAPG2(1:N) = nan;
% differences_objective_IPG(1:b) = nan;

ax = subplot(211);
hold off
semilogy([its_IAPG(2:end).f]);
hold on
plot(repmat(SVM_f(Bstar_cent), 1, max_its-1))
% semilogy(differences_objective_IPG, 'linewidth', 2);
% hold on
% plot(differences_objective_IAPG, 'linewidth', 2);
% plot(differences_objective_IAPG2, 'linewidth', 2);
grid on
ylabel('|f(x)-f(x*)|');
legend('IPG', 'IAPG', 'IAPG-light');
title('Objective');
stepsizes = [its_IAPG(1:end-1).x] - [its_IAPG(2:end).x];
stepsize_IAPG = arrayfun(@(i) norm(stepsizes(:, i)), 1:length(stepsizes));

ax2 = subplot(212);
linkaxes([ax, ax2], 'x');
hold off
semilogy(stepsize_IAPG);
grid on 
xlabel('Iterations');
ylabel('|| x_{k} - x_{k+1} ||');
title('Stepsize');


%% Plot feasibility percentage
initfig('Feasibility', 2);
ax = subplot(211);
% plot([its_IPG.feas]);
plot([its_IAPG.feas]);
grid on;
% plot([its_IAPG2.feas]);
% legend('IPG', 'IAPG', 'IAPG-light');
ylabel('% violated');

norm_subgrads = arrayfun(@(i) norm([its_IAPG(i).subgrad]), 1:max_its);

ax2 = subplot(212);
linkaxes([ax, ax2], 'x');
semilogy(norm_subgrads);
h = ylabel('$\| \tilde \nabla h_i(x_{k+1})\|$');
set(h, 'interpreter','latex');
xlabel('iteration');


%% plot times
initfig('Times', 3);
% plot([its_IPG.time]);
plot([its_IAPG.time]);
% plot([its_IAPG2.time]);
legend('IPG', 'IAPG', 'IAPG-light');
ylabel('Time per iteration');
xlabel('Iteration');