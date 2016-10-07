%% Image code for expressing sensitivity of m
% (c) OtH 13-5-16
% plots the number of samples needed for different eps, beta and d
%%

clf
betas = logspace(-2,-10, 20);            % confidence level
epss = linspace(0.1, 0.05, 20);          % violation level
ds = [10 50 100];                        % dimension of x
N = zeros(length(betas), length(eps), length(ds));

for i = 1:length(ds)
    d = ds(i);
    for j = 1:length(epss)
        eps = epss(j);
        N(:, j, i) = 2/eps*(d-1+log(1./betas));
    end
    surf(betas, epss, N(:,:,i)', 'DisplayName', sprintf('$\\zeta$ = %i', d));
    hold on
end
hold on
grid on
xlabel('\beta')
ylabel('\epsilon')
zlabel('N')
% h = legend('show', 'location', 'best')
% set(h, 'interpreter', 'latex')
set(gca, 'xscale', 'log');
% set(gca, 'yscale', 'log');
campos([0.0509109356937080,0.366832719972212,24779.6571355342]);
title('Effect of \epsilon, \beta and \zeta on required sample size N')
% xticklabels = get(gca, 'XTickLabels');
% set(gca, 'XTicklabels', ['10^{-10}', xticklabels(2:end)']')
%% vary d and beta
clf
betas = logspace(-2,-10, 20);            % confidence level
epss = [0.1 0.075 0.05];                 % violation level
ds = linspace(5, 500, 20);               % dimension of x
N = zeros(length(betas), length(ds), length(epss));
for i = 1:length(epss)
    eps = epss(i);
    for j = 1:length(ds)
        d = ds(j);
        N(:, j, i) = 2/eps*(d-1+log(1./betas));
    end
    surf(betas, ds, N(:,:,i)', 'DisplayName', sprintf('eps = %g', eps));
    hold on
end
hold on
grid on
xlabel('\beta')
ylabel('d')
zlabel('N')
legend('show', 'location', 'best')
set(gca, 'xscale', 'log');
% set(gca, 'yscale', 'log');
campos([0.0717461869828726,0.340461572134018,10019.1866590218]);
title('Effect of \epsilon, \beta and d on required sample size N')
% xticklabels = get(gca, 'XTickLabels');
% set(gca, 'XTicklabels', ['10^{-10}', xticklabels(2:end)']')
