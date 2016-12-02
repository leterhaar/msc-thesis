%% Initialize
addpath('../misc');
clear
yalmip('clear');
d = 2;
N = 50;
xs = randn(N, d);
ys = zeros(N,1);
ys(xs(:,1) < 0) = -1;
ys(xs(:,1) > 0) = 1;
b = 5;
%% Create problem

B = sdpvar(d,1, 'full');
Obj = 0.5*norm(B,2)^2;
C_cent = [];

for i = 1:N
    C_cent = [C_cent, ys(i) * xs(i, :) * B >= 1];
end

opt_settings = sdpsettings('verbose', 0);
diagnostics = optimize(C_cent, Obj, opt_settings);
assert(not(diagnostics.problem));
Bstar = value(B);

% define B0, somewhere in the feasible set
diagnostics = optimize(C_cent, [], opt_settings);
assert(not(diagnostics.problem));
B0 = value(B);

%% plot 
initfig('IAPG SVM', 1);
plot(xs(ys == -1, 1), xs(ys == -1, 2), 'x');
plot(xs(ys == 1, 1), xs(ys == 1, 2), '+');
the_axis = axis;
% plot centralized solution
xcoords = [-10 10];
plot(Bstar(2)*xcoords, Bstar(1)*xcoords, '--', 'linewidth', 2);
h = animatedline;
h.LineWidth = 2;
h.LineStyle = '--';
axis(the_axis)
addpoints(h, B0(2)*xcoords, B0(1)*xcoords);
point = animatedline;
point.Marker = '*';
point.MarkerEdgeColor = 'red';
point.MarkerSize = 15;
drawnow


%%
k = 0;
x0 = B0;
xk = zeros_like(x0);
z = zeros_like(x0);
all_i = 1:N;
max_its = 50;

x_sdp = sdpvar(d,1, 'full');
C_det = [];

its = struct(   'x', [], ...        value for x
                'f', [], ...        objective function value
                'subgrad', []...    subgradient at x
                ); 
                
opt_settings = sdpsettings('verbose', 0);


%% prepare first N iterations
[its(1:N).x] = deal(x0);
[its(1:N).f] = deal(0.5*norm(x0)^2);
[its(1:N).subgrad] = deal(zeros_like(x0));

%% start main iterations
assert(length(its) == N, 'Please prepare first');

p = progress('Iterating', max_its-N);
for k = N:max_its
    
    % select the next agent
    xk = its(k).x;
    ak = 1/(k+1);
    i = rem(k-1, N) + 1;
    
    
    past_gradients = zeros_like(x0);
    past_subgrads = zeros_like(x0);
    
    % loop over other indices than i
    for j = all_i(all_i ~= i)
        
        % for each i, find the previous k that it was active
        l = floor(k/N)*N - (N-j);
        
        past_gradients = past_gradients + its(l).x;
        past_subgrads = past_subgrads + its(l).subgrad;
    end
        
    z = xk - ak * (its(k).x + past_gradients);
    
    % project on feasible set i
    C_scen = ys(i) * xs(i,:) * x_sdp >= 1;
    Obj = past_subgrads' * (x_sdp - z) + 1/(2*ak) * norm(x_sdp - z, 2)^2;
    
    
    info = optimize([C_det, C_scen], Obj, opt_settings);
    assert(not(info.problem), sprintf('Problem optimizing: %s', info.info));

%    x2 = prox(Obj, ak, z, x_sdp, C_scen);
%     assert(all_close(x, x2));


    % store x(k+1) and its gradient and objective function
    x = value(x_sdp);
    its(k+1).x = x;
    its(k+1).f = 0.5*norm(x)^2;
    its(k+1).subgrad = x-z;
    
    % draw
%     clearpoints(h);
%     clearpoints(point);
%     addpoints(h, x(2)*xcoords, x(1)*xcoords);
%     plot(x(2)*xcoords, x(1)*xcoords, '--', 'Color', [0.8 0.8 0.8]);
%     addpoints(point, xs(i, 1), xs(i, 2));
%     drawnow
    
    p.ping(); 
    
end

%% plot second figure
initfig('IAPG SVM Objective', 2);
hold off
semilogy([its(N:end).f], 'linewidth', 2);
hold on
plot(repmat(0.5*norm(Bstar)^2, 1, max_its-N), 'linewidth', 2)
grid on
ylabel('|f(x)- f(x*)|');
xlabel('Iteration');

%%
initfig('IAPG SVM Solution', 3);
subplot(311);
xs_its = [its.x]';
hold on
grid on
plot(xs_its(N+1:end, 1), 'LineWidth', 2);
plot(repmat(Bstar(1), 1, max_its-N), 'linewidth', 2);

subplot(312);
hold on
grid on
plot(xs_its(N+1:end, 2), 'LineWidth', 2);
plot(repmat(Bstar(2), 1, max_its-N), 'linewidth', 2);

diff = xs_its(N+1:end, :) - repmat(Bstar', length(xs_its)-N, 1);
norms = arrayfun(@(x) norm(diff(x, :)), 1:length(xs_its)-N);
subplot(313);
semilogy(norms, 'linewidth', 2);
hold on
grid on

