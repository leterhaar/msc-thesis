%% quick and dirty plot to compare Js
figure(3);
set(gcf, 'name', 'ACC vs cACC');
assert(exist('Js_acc', 'var') == 1, 'No acc results available');
assert(exist('Js', 'var') == 1, 'No results available');
assert(exist('Obj', 'var') == 1, 'No centralized result available');
clf
Jc = value(Obj);

hold on
grid on
xlabel('Iteration');
ylabel('J(x*)');
title('Objective in time w/ and w/o consensus ACC');
h_without = plot(-1, Jc, '-', 'LineWidth', 2);
h_with = plot(-1, Jc, ':', 'LineWidth', 2);
h_central = plot(-1, Jc, '-.', 'linewidth', 1.2);
plot(1:t, value(Obj)*ones(1,t), 'linewidth', 1.2, 'color',h_central.Color,...
                                         'linestyle', h_central.LineStyle);
plot(Js, 'color', h_with.Color, 'LineWidth', 2, ...    
                                            'LineStyle', h_with.LineStyle);
plot(Js_acc, 'color', h_without.Color, 'LineWidth', 2,...
                                         'LineStyle', h_without.LineStyle);

set(gca, 'xtick', 1:t);
xlim([1 t]);
set(gca, 'xticklabels', (0:t-1))
legend('Agents (ACC)', 'Agents (cACC)', 'Centralized', 'location', 'se');
