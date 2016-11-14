load('../results/DC_scenarios_case30.mat');

m = length(agents);
T = agents{1}.k;

disagreements = zeros(T, m);
fieldnames = {'P_G', 'R_us', 'R_ds', 'd_us', 'd_ds'};
centralX = [];
for i = 1:5
    centralX = [centralX; central.X.(fieldnames{i})];
end

for i = 1:m
    
    % load agent
    agent = agents{i};
    
    % store Xs in one vector over time 
    Xs = [];
    for j = 1:5
        Xs = [Xs; agent.X.(fieldnames{j})];
    end
    
    % calculate difference between central solution and local iterates
    for t = 1:T
        disagreements(t, i) = norm(Xs(:, t) - centralX);
    end
    
end

%% plot
initfig('Disagreement CPA', 1);
hold off
semilogy(disagreements, 'linewidth', 1.5, 'color', blue);
xlabel('Iteration')
title('Disagreement per iteration of CP algorithm');
ylabel('|| x_i - z ||')
xlim([1 T]);
ylim([0.5 10]);
set(gca, 'YMinorTick', 'on');
grid on

