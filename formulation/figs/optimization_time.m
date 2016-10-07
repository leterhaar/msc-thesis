addpath('../../misc');
load('../results/results_P1P3s_mosek.mat');

networks = {'case14', 'case_ieee30'};
network_titles = {'14 bus', '30 bus'};
N = length(networks);
h = figure(1);
set(h, 'name', 'Time');
clf


for i = 1:N
    % select the right data
    results_network = results(strcmp({results.network}, networks{i}));
    n = length(results_network);
    results_network_P1 = results_network(strcmp({results_network.formulation}, 'P1'));
    results_network_P3 = results_network(strcmp({results_network.formulation}, 'P3'));
    results_network_P3s = results_network(strcmp({results_network.formulation}, 'P3*'));
    
    % set up subplot
    subplot(N,1,i);
    cla
    grid on
    hold on
    
    % plot
    h = plot([results_network_P1.scenarios], [results_network_P1.opt_time], '-s');
    h.MarkerFaceColor = h.Color;
    h = plot([results_network_P3.scenarios], [results_network_P3.opt_time], '-^');
    set(h, 'markerfacecolor', get(h, 'color'));
    h = plot([results_network_P3s.scenarios], [results_network_P3s.opt_time], '-p');
    set(h, 'markerfacecolor', get(h, 'color'));
    legend('P1', 'P3', 'P3*', 'location', 'nw');
    title(['Optimization time for ' network_titles{i} ' network']);
    ylabel('Time [s]');
    
    set(gca, 'YScale', 'log');

end


xlabel('Number of scenarios');
dock
