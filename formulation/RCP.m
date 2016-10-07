%% A Random Convex Program

% settings
n = 2;
N = 10;

% scenario generation
delta_1 = rand(n, N);
delta_2 = randn(n, N);
delta_3 = randn(1, N);

% define decision variable and objective function
x = sdpvar(n,1,'full');
obj = x' * x;

% define constraint
C = [];
figure
grid on; hold on;
for i = 1:N
    C = [C, x' * diag(delta_1(:, i)) * x + delta_2(:, i)' * x ...
                            + delta_3(i) >= x];
    for j = 1:n
        xplot = linspace(-10, 10, 1000);
        subplot(n,1,j)
        hold on
        plot(xplot, delta_1(j,i) .* xplot.^2 + delta_2(j, i) .* xplot ...
                                        + repmat(delta_3(i),1, 1000));
    end
end

% optimize
settings = sdpsettings('verbose', 1);
optimize(C, obj, settings);