load('../results/DC_scenarios_case30.mat');
addpath('../../misc');

% prepare fig
figure(5);
set(gcf, 'name', 'Consensus');
clf
dock
grid on
hold on


N = length(agents{1}.X)-1;                % number of iterations
Na = length(agents);                    % number of agents
Ng = length(agents{1}.X(1).P_G);        % number of generators


% plot for every agent the first, halfway and end dispatch
P_Gs = zeros(Ng, Na, 2);
R_us = zeros(Ng, Na, 2);
R_ds = zeros(Ng, Na, 2);

% loop over agents
for a = 1:Na
    agent = agents{a};
    
    % store P_G, R_us and R_ds
    P_Gs(:,a,1) = agent.X(1).P_G;
    P_Gs(:,a,2) = agent.X(N).P_G;
    R_us(:,a,1) = agent.X(1).R_us;   
    R_us(:,a,2) = agent.X(N).R_us;
    R_ds(:,a,1) = -agent.X(1).R_ds;
    R_ds(:,a,2) = -agent.X(N).R_ds;
end

titles = {'Generator dispatch and reserve requirements after 1 iteration', sprintf('Generator dispatch and reserve requirements after %i iterations', N)};

subplot(311);
hold on
bar(-1, 0, 'b');
bar(-2, 0, 'r');
bar(-3, 0, 'g');
bar(-4, 0, 'b', 'facealpha', 0.3, 'linestyle', ':');
bar(-4, 0, 'r', 'facealpha', 0.3, 'linestyle', ':');
bar(-4, 0, 'g', 'facealpha', 0.3, 'linestyle', ':');

legend('P_G per agent', 'R_{us} per agent', 'R_{ds} per agent', ...
    'P_G centralized', 'R_{us} centralized', 'R_{ds} centralized',  ...
    'location', 'neo');

for i = 1:2
    subplot(2,1,i);
    title(titles{i});
    grid on
    hold on
    bar(central.X.P_G, 'b', 'facealpha', 0.3, 'linestyle', ':');
    bar(central.X.R_us, 'r', 'facealpha', 0.3, 'linestyle', ':');
    bar(-central.X.R_ds, 'g', 'facealpha', 0.3, 'linestyle', ':');
    bar(P_Gs(:,:,i), 'b');
    bar(R_us(:,:,i), 'r');
    bar(R_ds(:,:,i), 'g');    
    ylabel('Power');
    set(gca, 'xtick', 1:Ng);
    set(gca, 'xlim', [0.5 Ng+0.5]);
end
    xlabel('Generator');
    
%% iterations on x axis

figure(6);
set(gcf, 'name', 'Iterations');
dock;
clf;

N_k = agents{1}.k-1;
N_G = size(agents{1}.X(1).P_G, 1);

ylabels = {'P_G', 'R_{us}' , 'R_{ds}'};
keys = {'P_G', 'R_us', 'R_ds'};
titles = {'Generator dispatch' , 'Upspinning reserve', 'Downspinning reserve'};
for i = 1:3
    
    subplot(3,1,i);
    hold on;
    grid on;
    title(titles{i});
    ylabel(ylabels{i});
    
    for a = 1:Na
        set(gca, 'ColorOrderIndex', 1);
        plot([agents{a}.X.(keys{i})]', '-.');
    end
    
    % plot central solution    
    set(gca, 'ColorOrderIndex', 1);
    plot(repmat(central.X.(keys{i}), 1, N_k)', '-');
    
    xlim([1, N_k]);

    
end
xlabel('Iteration');

