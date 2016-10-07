%% show upper bound for Helly dimension for different definitions

% assume N_l = 4/3 * N_b
% assume N_t = 24

N_b = linspace(0, 350, 1000);
N_t = 24;
N_l = 4/3*N_b;
epsilon = 0.1;
beta = 1e-5;

zeta_n = (10*N_b+4*N_b.^2)*N_t;
zeta_g = (4*N_b+3*N_l+4)*N_t*(2*N_t+1);
N_n = ceil(2/epsilon*(zeta_n-1+log(1/beta)));
N_g = ceil(2/epsilon*(zeta_g-1+log(1/beta)));

clf
dock
grid on
hold on
plot(N_b, zeta_n, '-.', 'linewidth', 3, 'DisplayName', 'd');
plot(N_b, zeta_g, '-', 'linewidth', 3, 'DisplayName', 'r(m+1)');
legend('show', 'location', 'northwest');
xlabel('Number of buses N_b');
ylabel('Upper bound on Helly-dimension \zeta');
title('Relation between network size and bounds on \zeta and N');
ylim([0 max([zeta_n, zeta_g])]);

yyaxis right
% plot(N_b, N_n, '-.', 'linewidth', 1, 'DisplayName', 'd');
% plot(N_b, N_g, '-', 'linewidth', 1, 'DisplayName', 'r(m+1)');
ylim([0 max([N_n, N_g])]);
set(gca, 'YColor', [0 0 0]);
ylabel('Number of scenarios required N');
%
crossing_idx = find(zeta_n > zeta_g, 1);
x = N_b(crossing_idx);
y = zeta_n(crossing_idx);
% prettyText(x,y+5e5, sprintf('d > r(m+1) for N_b > %.f', x), 'northwest');