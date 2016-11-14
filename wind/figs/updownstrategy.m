initfig('UD strat', 1);

us = [0.15 0.5 0.35];
ds = [0.15 0.1 0.75];

as = area([0 1],[zeros(1,3); -us]);
colors = [green; blue; brown];
opacities = [0.8 0.8 0.8];
for i = 1:3
    a = as(i);
    a.FaceColor = colors(i, :);
%     a.FaceAlpha = opacities(i);
end
as = area([-1 0],[ds; zeros(1,3)]);
reverse = [3 2 1];
for i = 1:3
    a = as(i);
    a.FaceColor = colors(reverse(i), :);
%     a.FaceAlpha = opacities(i);
end

legend('G1', 'G2', 'G3', 'location', 'ne');
set(gca, 'XTickLabels', [], 'YTickLabels', []);