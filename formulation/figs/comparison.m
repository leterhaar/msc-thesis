% load results
addpath('../../misc');

% P1
load('../results/results_P1.mat')
P1_opt = opt_times;
P1_total = total_times;

% P2
load('../results/results_P3.mat')
P2_opt = opt_times;
P2_total = total_times;

% P2star
load('../results/results_P3star.mat')
P2star_opt = opt_times;
P2star_total = total_times;

figh = initfig('Comparison', 1);
% hold off
% plot(network_sizes, P1_opt, '+--', 'color', blue,  'markerface', blue, 'linewidth', 2)
% hold on
% grid on
hold off
loglog(network_sizes, P1_total, 'o-', 'color', blue, 'markerface', blue, 'linewidth', 2)
grid on
hold on
% plot(network_sizes, P2_opt, '+--', 'color', green, 'markerface', green,  'linewidth', 2)
plot(network_sizes, P2_total, 'o-', 'color', green,  'markerface', green, 'linewidth', 2)
% plot(network_sizes, P2star_opt, '+--', 'color', orange,  'markerface', orange, 'linewidth', 2)
plot(network_sizes, P2star_total, 'o-', 'color', orange,  'markerface', orange, 'linewidth', 2)
% xlim([network_sizes(1) network_sizes(end)]);

X = [ones(length(network_sizes),1) network_sizes'];
b = X\P2_total'
t = [5 300 5000 24*5000];
plot(t, b(2)*t+b(1), 's-r');