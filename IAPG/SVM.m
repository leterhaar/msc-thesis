%% Initialize
addpath('../misc');
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

B_opt = value(B);
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
max_its = 500;

if max_its > N
    random_i = [];
    for j = 1:ceil(max_its/N)
        random_i = [random_i; randperm(N, N)];
    end
    random_i = random_i(1:max_its);
else
    random_i = randperm(N, max_its);
end

x_sdp = sdpvar(d,1, 'full');
C_det = [];

its = struct(   'x', [], ...        value for x
                'f', [], ...        objective function value
                'grad', []...      gradient at x
                ); ...   general infeasibility
                    
its(1).x = x0;
its(1).grad = x0;
its(1).f = 0.5*norm(x0)^2;
% its(1).check = DC_check(x0, dc, wind);

opt_settings = sdpsettings('verbose', 0);

%% prepare first b iterations
p = progress('Preparing', b);
for k = 1:b
    xk = its(k).x;
    ak = 1/(k+1);
    i = random_i(k);
    
    past_gradients = zeros_like(x0);
    
    % loop over other indices than i
    for j = all_i(all_i ~= i)
        if k > 1
            l = randi(k-1); % let the delay be dependent (for preparation only)
        else
            l = 0;
        end
        xl = its(k-l).x;
        past_gradients = past_gradients + xl;
    end
        
    z = xk - ak * (its(k).x + past_gradients);
    
    % project on feasible set i
    C_scen = ys(i) * xs(i,:) * x_sdp >= 1;
    Obj = 1/(2*ak) * norm(x_sdp - z, 2)^2;
    
    info = optimize([C_det, C_scen], Obj, opt_settings);
    assert(not(info.problem), sprintf('Problem optimizing: %s', info.info));
    
    % store x(k+1) and its gradient and objective function
    x = value(x_sdp);
    its(k+1).x = x;
    its(k+1).f = 0.5*norm(x)^2;
    its(k+1).grad = x;
    
    % draw
    clearpoints(h);
    clearpoints(point);
    addpoints(h, x(2)*xcoords, x(1)*xcoords);
    plot(x(2)*xcoords, x(1)*xcoords, '--', 'Color', [0.8 0.8 0.8]);
    addpoints(point, xs(i, 1), xs(i, 2));
    drawnow
    
    p.ping(); 
end

%% start main iterations
assert(length(its) > b, 'Please prepare first');

p = progress('Iterating', max_its-b);
try
for k = b+1:max_its
    
    xk = its(k).x;
    ak = 1/(k+1);
    i = random_i(k);
    
    past_gradients = zeros_like(x0);
    
    % loop over other indices than i
    for j = all_i(all_i ~= i)
        % choose a random delay
        l = randi(b);
        past_gradients = past_gradients + its(k-l).x;
    end
        
    z = xk - ak * (its(k).x + past_gradients);
    
    % project on feasible set i
    C_scen = ys(i) * xs(i,:) * x_sdp >= 1;
    Obj = 1/(2*ak) * norm(x_sdp - z, 2)^2;
    
    info = optimize([C_det, C_scen], Obj, opt_settings);
    assert(not(info.problem), sprintf('Problem optimizing: %s', info.info));
    
    % store x(k+1) and its gradient and objective function
    x = value(x_sdp);
    its(k+1).x = x;
    its(k+1).f = 0.5*norm(x)^2;
    its(k+1).grad = x;
    
    % draw
%     clearpoints(h);
%     clearpoints(point);
%     addpoints(h, x(2)*xcoords, x(1)*xcoords);
%     plot(x(2)*xcoords, x(1)*xcoords, '--', 'Color', [0.8 0.8 0.8]);
%     addpoints(point, xs(i, 1), xs(i, 2));
%     drawnow
    
    p.ping(); 
    
end
catch me
    warning(me.getReport())
end


%% plot second figure
initfig('IPG SVM Objective', 4);
hold off
semilogy([its(N:end).f], 'linewidth', 2);
hold on
plot(repmat(0.5*norm(Bstar)^2, 1, max_its-N), 'linewidth', 2)
grid on
ylabel('|f(x)- f(x*)|');
xlabel('Iteration');

%%
initfig('IPG SVM Solution', 5);
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

