t = 0:0.1:24;
y = 1+-0.2*sin(t/12*pi);
clf
hold on
grid on
plot(t, ones(length(t),1), 'k-.', 'LineWidth', 2, 'color', 0.6*ones(3,1));
plot(t, y, 'k', 'linewidth', 2);
area(t(y<1), y(y<1), 1, 'FaceColor', [216,179,101]./255, 'EdgeColor', 'none')
area(t(y>1), y(y>1), 1, 'FaceColor', [90,180,172]./255, 'EdgeColor', 'none')

plot(t, ones(length(t),1), 'k-.', 'LineWidth', 2, 'color', 0.6*ones(3,1));
plot(t, y, 'k', 'linewidth', 2);
legend('Scheduled production', 'Actual consumption',  ...
       'Downspinning required', 'Upspinning required');

set(gca, 'YTickLabel', {});
set(gca, 'XTickLabel', {});

xlabel('Time');
ylabel('Power');
title('Illustration of up- and downspinning regions')
ylim([0.5, 1.5]);
xlim([0,24]);
