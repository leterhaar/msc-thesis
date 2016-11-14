addpath('../../misc/');

h = initfig('PF',1);

set(gca, 'FontName', 'Fira Sans');

t = 0:0.1:10*pi;
u = cos(t);
plot(t, u, 'LineWidth', 2, 'color', blue);
xlim([0 10*pi]);