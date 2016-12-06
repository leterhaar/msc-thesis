%% plot objectives
addpath('../../misc');

initfig('Objectives', 1);
differences_objective_IPG = abs([its_IPG(2:end).f] - SVM_f(Bstar_cent));
differences_objective_IAPG = abs([its_IAPG(2:end).f] - SVM_f(Bstar_cent));
differences_objective_IAPG2 = abs([its_IAPG2(2:end).f] - SVM_f(Bstar_cent));

ax = subplot(211);
hold off
semilogy(differences_objective_IPG, 'linewidth', 2);
hold on
plot(differences_objective_IAPG, 'linewidth', 2);
plot(differences_objective_IAPG2, 'linewidth', 2);
grid on
ylabel('|f(x)-f(x*)|');
legend('IPG', 'IAPG', 'IAPG-light');
title('Objective');


% Plot feasibility percentage
ax2 = subplot(212);
plot([its_IPG.feas]);
hold on
title('Feasibility');

plot([its_IAPG.feas]);
grid on;
plot([its_IAPG2.feas]);
legend('IPG', 'IAPG', 'IAPG-light');
ylabel('% violated');
linkaxes([ax, ax2], 'x');
%% plot times
initfig('Times', 2);
plot([its_IPG.time]);
plot([its_IAPG.time]);
plot([its_IAPG2.time]);
legend('IPG', 'IAPG', 'IAPG-light');
ylabel('Time per iteration');
xlabel('Iteration');