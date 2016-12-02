%% initialize models
% add paths
addpath('../networks');
addpath('../misc');
addpath('../formulation_dc');
addpath('../wind');
yalmip('clear');

dc = DC_model('case14a');
dc.set_WPG_bus(9);

N_t =   24;     % number of timesteps
N =     200;    % number of scenarios
b = 5;          % maximum delay
max_its = 100;  % maximum number of iterations

wind = wind_model(dc, N_t, 0.2);
wind.generate(N);
%%
k = 0;
x0 = zeros(5*dc.N_G, N_t);
xk = zeros_like(x0);
z = zeros_like(x0);
all_i = 1:N;

if max_its > N
    random_i = [];
    for j = 1:ceil(max_its/N)
        random_i = [random_i; randperm(N, N)];
    end
    random_i = random_i(1:max_its);
else
    random_i = randperm(N, max_its);
end

x_sdp = sdpvar(5*dc.N_G, N_t, 'full');
C_det = DC_cons_det(x_sdp, dc, wind);
C_scen = [];
for i = 1:N
    C_scen = [C_scen, DC_cons_scen(x_sdp, dc, wind.slice(i))];
end
opt_settings = sdpsettings('verbose', 0);

diagnostics = optimize([C_det, C_scen], [], opt_settings);
assert(not(diagnostics.problem), 'Problem optimizing');
x0 = value(x_sdp);

its = struct(   'x', [], ...        value for x
                'f', [], ...        objective function value
                'grad', [],...      gradient at x
                'check', []); ...   general infeasibility
                    
its(1).x = x0;
its(1).grad = DC_gradient_f(x0, dc, wind);
its(1).f = DC_f(x0, dc, wind);
its(1).check = DC_check(x0, dc, wind);


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
        past_gradients = past_gradients + DC_gradient_f(xl, dc, wind);
    end
        
    z = xk - ak * (its(k).grad + past_gradients);
    
    % project on feasible set i
    C_scen = DC_cons_scen(x_sdp, dc, wind.slice(i));
    Obj = 1/(2*ak) * norm(x_sdp - z, 'fro')^2;
    
    info = optimize([C_det, C_scen], Obj, opt_settings);
    assert(not(info.problem), sprintf('Problem optimizing: %s', info.info));
    
    % store x(k+1) and its gradient and objective function
    its(k+1).x = value(x_sdp);
    its(k+1).f = DC_f(its(k+1).x, dc, wind);
    its(k+1).check = DC_check(its(k+1).x, dc, wind);
    its(k+1).grad = DC_gradient_f(its(k+1).x, dc, wind);
    p.ping(); 
end

%% start main iterations
assert(length(its) > b, 'Please prepare first');

p = progress('Iterating', max_its-b);
for k = b+1:max_its
    
    xk = its(k).x;
    ak = 1/(k+1);
    i = random_i(k);
    
    past_gradients = zeros_like(x0);
    
    % loop over other indices than i
    for j = all_i(all_i ~= i)
        % choose a random delay
        l = randi(b);
        past_gradients = past_gradients + its(k-l).grad;
    end
        
    z = xk - ak * (its(k).grad + past_gradients);
    
    % project on feasible set i
    C_scen = DC_cons_scen(x_sdp, dc, wind.slice(i));
    Obj = 1/(2*ak) * norm(x_sdp - z, 'fro')^2;
    
    info = optimize([C_det, C_scen], Obj, opt_settings);
    assert(not(info.problem), sprintf('Problem optimizing: %s', info.info));
    
    % store x(k+1) and its gradient and objective function
    its(k+1).x = value(x_sdp);
    its(k+1).f = DC_f(its(k+1).x, dc, wind);
    its(k+1).check = DC_check(its(k+1).x, dc, wind);
    its(k+1).grad = DC_gradient_f(its(k+1).x, dc, wind);
    p.ping(); 
    
end

%% calculate centralized solution
p = progress('Central', 2);
C_all = DC_cons_scen(x_sdp, dc, wind);
Obj = DC_f(x_sdp, dc, wind);
p.ping();
info = optimize([C_det, C_all], Obj, opt_settings);

assert(not(info.problem), sprintf('Problem optimizing: %s', info.info));
xstar = value(x_sdp);
fstar = DC_f(xstar, dc, wind);
p.ping();

%% plot solutions
plot_its = min(max_its, 250);
initfig('Solution', 1);
% xs = [its.x];
diffs = nan(plot_its);

for i = 1:plot_its
    diffs(i) = norm(its(i).x-xstar, 'fro');
end

semilogy(diffs, 'linewidth', 2);
hold on;
grid on;
xlabel('Iteration');
ylabel('|x-x*|');
title('(Lack of) Convergence of IAPG algorithm');
ylims = ylim;
ylim([0 ylims(2)]);
%% plot objective value

initfig('Objective', 2);
fs = abs([its(2:plot_its).f] - fstar);

semilogy(fs, 'linewidth', 2);
hold on
grid on
xlabel('Iteration')
ylabel('|f(x) - f(x*)|');
title('Objective values of IAPG algorithm');


ylims = ylim;
ylim([0 ylims(2)]);
