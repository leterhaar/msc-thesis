%% Plot wind scenarios used for optimization and validation
figure(3)
clf
dock
hold on
grid on

h = area(1:24, [min(P_wscenssim, [], 2) max(P_wscenssim, [], 2)]);
h(1).FaceColor = 'none';
h(1).EdgeColor = 'none';
h(2).EdgeColor = 'none';
h(2).FaceColor = [199,234,229]./255;
plot(1:24, P_wscens, '-', 'Color', [1,133,113]./255, 'LineWidth', 0.8);
plot(1:24, P_wf, '-', 'Color', [166,97,26]./255, 'LineWidth', 2.4);
legend('Optimization', 'Validation', ...
        'Forecasted');
% NB fix legend labels in illustrator
set(gca, 'xtick', 1:24);
xlim([1,24]);
xlabel('Time [h]');
ylabel('Wind-power [MW]');
title('Wind-power trajectories over time');