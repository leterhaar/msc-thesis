clf
hold on

% which hour to plot
h = 2;

% results
bars = [P_Goptdet(:, h) P_Gopt(:, h) R_us(:, h) R_ds(:, h)];

hold on
grid on
bar(bars);
xlabel('Generator index');
ylabel('P^G [MW]');
title(sprintf('DC OPF at hour %d, forecasted wind = %.2g MWh', h, wind.P_wf(h)));

legend('P^G (without wind)', 'P^G (with wind)', 'R^{us}', 'R^{ds}');
set(gca, 'xticklabels', {'', 1, 2, 3, 4, 5, 6, ''});