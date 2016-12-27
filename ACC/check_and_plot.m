function [f1, f2, f3] = check_and_plot(agents, optimal_objective, residual_fcn)
%% [f1, f2, f3] = check_and_plot(agents, optimal_objective, C_all)
% checks feasibity and plots convergence for iterations of ACC(A) algorithm
% 
% returns figure handles for adding titles etc

    %% calculate convergence and feasibility
    K = length(agents(1).iterations);
    N = length(agents);
    convergence = nan(K,N);
    feasibility = nan(K,N);
    time_per_iteration = nan(K,N);
    optimizations_run = zeros(K,1);
    no_cons_used = nan(K,N);
    timing = zeros(K,1);

    p = progress('Checking constraints', N);

    for i = 1:N
        for k = 1:K
            % calculate difference with centralized objective
            convergence(k,i) = abs(agents(i).iterations(k).J ...
                                                    - optimal_objective);

            % calculate feasibility percentage
            residuals = residual_fcn(agents(i).iterations(k).x);
            feasibility(k,i) = sum(residuals < -1e-6) / length(residuals) * 100;

            % store times
            time_per_iteration(k,i) = agents(i).iterations(k).time;
            timing(k) = timing(k) + agents(i).iterations(k).time;


            % store total number of iterations run
            if k > 1
                optimizations_run(k) = optimizations_run(k) + ...
                                        agents(i).iterations(k).info.optimized;

                no_cons_used(k,i) = agents(i).iterations(k).info.num_cons;

            end
        end
        p.ping();
    end
    
    verify(sum(feasibility(end,:)) == 0, 'Not all feasible');

    %% plot
    f1 = initfig('ACC iterations', 4);
    ax = subplot(211, 'YScale', 'log');
    grid on
    hold on
    plot(convergence, 'color', blue);
    ylabel('|f(x_k^i) - f(x^*) |')

    ax2 = subplot(212);
    linkaxes([ax ax2], 'x');
    grid on
    hold on
    plot(feasibility, 'color', green);

    ylabel('% violated');
    xlabel('iterations');

    f2 = initfig('Timing', 2);
    ax = subplot(211);
    grid on
    hold on
    plot(time_per_iteration, 'o', 'color', green);
    plot(timing, 'o-', 'color', green);
    ylabel('Time per iteration');

    ax2 = subplot(212);
    linkaxes([ax ax2], 'x');
    grid on
    hold on
    plot(optimizations_run, 'color', green);
    ylabel('Optimizations run');
    xlabel('Iteration');

    f3 = initfig('No of constraints used', 3);
    hs1 = plot(no_cons_used, 'color', green);
    h2 = plot(repmat(length(residuals), K, 1), '--');
    legend([hs1(1) h2], '|A_i^{(k)}|', '|C_{all}|');
    uistack(f1);
    figure(f1);
end