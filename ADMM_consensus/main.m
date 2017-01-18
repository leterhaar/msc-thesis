clc
clear
yalmip('clear');

if not(exist('AC_model', 'file'))
    addpath('../wind');
    addpath('../misc');
    addpath('../networks');
    addpath('../formulation');
end

%% load models
N_t = 24;   % optimization horizon
N = 5;      % number of scenarios used for optimization
t = 1; % timestep used for this demonstration (todo: add for everything)
mu = 10;
max_its = 500;

% load network and wind models
ac = AC_model('case5');
ac.set_WPG_bus(2);
wind = wind_model(ac, N_t, 0.2);

% generate a number of scenarios
wind.dummy(N);
wind.use_extremes(t);
N = 2;
%% initialize agents
for i = 1:N
    agents(i) = ADMM_consensus_agent(ac, wind.slice(i), t, mu);
end

%% start iterating
maintimer = tic;
fprintf(' it# |    residual | time (s)\n');
fprintf('=============================\n');
residuals = nan(N,1);
averages = zeros(4*ac.N_b^2+4*ac.N_G,1);
k = 1;
while k <= max_its
    previous_averages = averages;
    averages = zeros(4*ac.N_b^2+4*ac.N_G,1);
    for i = 1:N
        residuals(i) = norm(previous_averages - agents(i).it(k+1).x)^2;
        averages = averages + agents(i).update(previous_averages);
    end    
    
    % calculate average
    averages = averages./N;

    residual = sum(residuals);
    fprintf('%4i | %.5e | %3.2f\n', k, residual, toc(maintimer));
    k = k + 1;
    
    if residual < 1e-4
        break;
    end
end

%% plot
initfig('ADMM iterations', 3);
hold off;
for i = 1:N
    semilogy([agents(i).it(:).residual]);
    hold on;
end
grid on; box on;
xlabel('iterations');
ylabel('Primal residuals $\| x_i^{(k)} - \bar x^{(k)} \|$', 'interpreter', 'latex');
%% extraxt xs and check feasibility
for i = 1:N
    fprintf('Checking agent %i\n', i);
    fprintf('=================\n');
    
    Wf = mat(agents(i).it(k).x(1:4*ac.N_b^2));
    dus = agents(i).it(k).x(4*ac.N_b^2+1:4*ac.N_b^2+ac.N_G);
    dds = agents(i).it(k).x(4*ac.N_b^2+ac.N_G+1:4*ac.N_b^2+2*ac.N_G);
    fprintf('Wf is_psd \t: %i\nWf rank\t\t: %i\n', is_psd(Wf), svd_rank(Wf, 10));
    
    [flag, msg] = simulate_network(ac, wind.slice(1), Wf, dus, dds, t, 1e-5);
    fprintf('Checked\t: %s\n\n', msg);
end