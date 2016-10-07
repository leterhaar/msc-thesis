rek = 2.5;
Ns = 500;
ss = linspace(rek*min(wind.P_m(t,:)), rek*max(wind.P_m(t,:)), Ns);

gammas = zeros(Ns, 2);
for i = 1:Ns
    C = Wf_opt + ss(i)*Wm_opt;
    gamma = eig(C);
    gammas(i,1) = gamma(1); % store the smallest
    gammas(i,2) = gamma(end); % store the biggest
end

pmin = ss(1)/rek;
pmax = ss(end)/rek;

clf
hold on
grid on
line = plot(ss, gammas(:,1), '-', 'LineWidth', 2);
% xlim([ss(1) ss(end)]);
ylabel('Smallest eigenvalue of W^s = W^f + P^m_i W^m');
xlabel('Wind power mismatch P^m');
xlims = [ss(1) ss(end)];
xlim(xlims);
ylims = ylim();
ylims = [ylims(1) ylims(2)+0.02];
dotted1 = plot([pmin pmin], [ylims(1)-1 ylims(2)+1], 'k--');
dotted2 = plot([pmax pmax], [ylims(1)-1 ylims(2)+1], 'k--');
xspacing = 0.1*(xlims(2)-xlims(1));
yspacing = 0.05*(ylims(2)-ylims(1));

not_psd = area([xlims(1) xlims(end)], ...
                [ylims(1)-1 ylims(1)-1]);
not_psd.FaceColor = 'r';
not_psd.FaceAlpha = 0.2;
not_psd.LineStyle = 'none';
psd = area([xlims(1)-1 xlims(2)+1], ...
                [ylims(2) ylims(2)+1]);
psd.FaceColor = 'g';
psd.FaceAlpha = 0.2;
psd.LineStyle = '-';
psd.EdgeColor = psd.FaceColor;
psd.LineWidth = 3;
psd.EdgeAlpha = 1;
axis([xlims ylims]);

uistack(dotted1, 'top');
uistack(dotted2, 'top');
uistack(line, 'top');

l = legend([psd not_psd], '$W^s \succeq 0$', '$W^s \not \succeq 0$');
l.Interpreter = 'latex';
l.FontSize = 12;

ax = gca;
xlabels = ax.XTickLabel;
xtick = ax.XTick;
xtick = [xtick pmin pmax];
xlabels{end+1} = 'P^m_{min}';
xlabels{end+1} = 'P^m_{max}';

% sort 
[xtick, order] = sort(xtick);
ax.XTick = xtick;
ax.XTickLabel = xlabels(order);

title('PSDness of W^s for different P^m');
