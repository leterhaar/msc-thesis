addpath('../../misc');
load('../results/results_P1P3_mosek.mat');

networks = {'case9','case14', 'case_ieee30'};
network_titles = {'9 bus', '14 bus', '30 bus'};
N = length(networks);
h = figure(3);
set(h, 'name', 'Rank');
clf


for i = 1:N
    % select the right data
    results_network = results(strcmp({results.network}, networks{i}));
    n = length(results_network);
    results_network_P1 = results_network(1:n/2);
    results_network_P3 = results_network(n/2+1:end);
    
    % set up subplot
    subplot(N,1,i);
    cla
    grid on
    hold on
    
    % plot
    h = plot([results_network_P1.scenarios], [results_network_P1.rank_perc], 's');
    set(h, 'MarkerFaceColor', get(h, 'color'));
    h = plot([results_network_P3.scenarios], [results_network_P3.rank_perc], '^');
    set(h, 'markerfacecolor', get(h, 'color'));
    legend('P1', 'P3', 'location', 'ne');
    title(['Rank 1 matrices for ' network_titles{i} ' network']);
    ylabel('Percentage');
    

end
xlabel('Number of scenarios');
dock