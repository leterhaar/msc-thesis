addpath('../../formulation');
load('../results/AC_dummy.mat');

% prepare fig
f = figure(5);
set(f, 'name', 'Consensus');
clf
grid on
hold on
dock

N = min([agents.k]);                % number of iterations
Na = length(agents);                    % number of agents


% plot for every agent the first, halfway and end dispatch
P_Gs = zeros(ac.N_G, Na, 3);
R_us = zeros(ac.N_G, Na, 3);
R_ds = zeros(ac.N_G, Na, 3);

% loop over agents
for a = 1:Na
    agent = agents(a);
    
    % store P_G
    for j = 1:ac.N_G
        k = ac.Gens(j);        
        P_Gs(j,a,1) = trace(ac.Y_k(k)*agent.X(1).W_f)+ac.P_D(k);
        P_Gs(j,a,2) = trace(ac.Y_k(k)*agent.X(N).W_f)+ac.P_D(k);
    end
    
    R_us(:,a,1) = agent.X(1).R_us;
    R_us(:,a,2) = agent.X(N).R_us;
    R_ds(:,a,1) = -agent.X(1).R_ds;
    R_ds(:,a,2) = -agent.X(N).R_ds;
end

% store P_Gs_central
P_G_central = zeros(ac.N_G, 1);
R_us_central = zeros(ac.N_G, 1);
R_ds_central = zeros(ac.N_G, 1);
for j = 1:ac.N_G
    k = ac.Gens(j);
    P_G_central(j) = trace(ac.Y_k(k)*central.X.W_f)+ac.P_D(k);
end

R_us_central = central.X.R_us;
R_ds_central = -central.X.R_ds;
titles = {'Generator dispatch and reserve requirements after 1 iteration', ...
    sprintf('Generator dispatch and reserve requirements after %i iterations', N)};


subplot(211);
hold on
bar(-1, 0, 'b');
bar(-2, 0, 'r');
bar(-3, 0, 'g');
bar(-4, 0, 'b', 'facealpha', 0.3, 'linestyle', ':');
bar(-4, 0, 'r', 'facealpha', 0.3, 'linestyle', ':');
bar(-4, 0, 'g', 'facealpha', 0.3, 'linestyle', ':');

legend('P_G per agent', 'R_{us} per agent', 'R_{ds} per agent', ...
    'P_G centralized', 'R_{us} centralized', 'R_{ds} centralized',  ...
    'location', 'nw');

for i = 1:2
    subplot(2,1,i);
    title(titles{i});
    grid on
    hold on
    bar(P_G_central, 'b', 'facealpha', 0.3, 'linestyle', ':');
    bar(R_us_central, 'r', 'facealpha', 0.3, 'linestyle', ':');
    bar(R_ds_central, 'g', 'facealpha', 0.3, 'linestyle', ':');
    bar(P_Gs(:,:,i), 'b');
    bar(R_us(:,:,i), 'r');
    bar(R_ds(:,:,i), 'g');    
    ylabel('Power');
    set(gca, 'xtick', 1:ac.N_G);
    set(gca, 'xlim', [0.5 ac.N_G+0.5]);
end
    xlabel('Generator');

%% iterations on x axis

figure(6);
set(gcf, 'name', 'Iterations');
dock;
clf;

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
        if strcmp(keys{i}, 'P_G')
            P_G = zeros(ac.N_G, N);
            for h = 1:N
                for j = 1:ac.N_G
                    k = ac.Gens(j);
                    P_G(j,h) = trace(ac.Y_k(k) * agents(a).X(h).W_f)+ac.P_D(k);
                end
            end
            plot(P_G', '-.');
        else
            plot([agents(a).X.(keys{i})]', '-.'); 
        end
    end
    
    % plot central solution    
    set(gca, 'ColorOrderIndex', 1);
    if strcmp(keys{i}, 'P_G')
            P_G = zeros(ac.N_G, N);
            for j = 1:ac.N_G
                k = ac.Gens(j);
                P_G(j,:) = trace(ac.Y_k(k) * central.X.W_f)+ac.P_D(k);
            end
            plot(P_G', '-', 'LineWidth', 1);
        else
            plot(repmat([central.X.(keys{i})], 1, N)', '-'); 
    end
    
    xlim([1, N]);

    
end
xlabel('Iteration');

