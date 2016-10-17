function format_result(m, wind, t, decided_vars, R)
% Formats the result of the AC OPF
%
% Arguments
% model : model used
% wind : wind model used
% t : current timestep
% decided_vars : cell with optimal values for {Wf, Ws, Rus, Rds, dus, dds}
% R : [N_G x N] matrix with reserve power

    if isa(m, 'AC_model')
        Wf_opt = decided_vars{1};
        Ws_opt = decided_vars{2};
        Rus_opt = decided_vars{3};
        Rds_opt = decided_vars{4};
        dus_opt = decided_vars{5};
        dds_opt = decided_vars{6};
        P_G = zeros(m.N_G,1);
    else
        P_G = decided_vars{1};
        Rus_opt = decided_vars{2};
        Rds_opt = decided_vars{3};
        dus_opt = decided_vars{4};
        dds_opt = decided_vars{5};
    end

    [~,N] = size(wind.P_w);
    % pick at much 100 samples to show
    N_samples = min(N, 15);
    format_samples = sort(randsample(N, N_samples))';
    



    % output table

    fprintf('\tGenBus \t\tP_G \td_us\td_ds \t R_us \t\t R_ds \t\t');
    for i = format_samples
        fprintf(' R_S%3i\t\t',i);
    end
    fprintf(['\n-------------------------------------------------------------' repmat('------------', 1, N_samples) '\n']);
    for k = 1:m.N_G
        busid = m.Gens(k);
        if isa(m, 'AC_model')
            P_G(k) = trace(m.Y_k(busid)*Wf_opt) + m.P_D(t,busid) - m.C_w(busid)*wind.P_wf(t);
        end
        fprintf(' %4i \t\t%f \t%.2f\t%.2f\t%f\t%f\t', busid, P_G(k), dus_opt(k), dds_opt(k), Rus_opt(k), Rds_opt(k));
        for i = format_samples
            fprintf('%f\t', R(k,i));
        end
        fprintf('\n');
    end
    fprintf(['-------------------------------------------------------------' repmat('------------', 1, N_samples) ' + \n']);
    fprintf('  TOTAL \t%f\t\t\t\t\t%f\t%f\t', sum(P_G), sum(Rus_opt), sum(Rds_opt));
    for i = format_samples
        fprintf('%f\t', sum(R(:,i)));
    end
    fprintf('\n\nP_D-P_wf\t%f\n  -P_m\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t', sum(m.P_D(t,:))-wind.P_wf(t));
    for i = format_samples
        fprintf('%f\t', -wind.P_m(t, i));
    end
    fprintf(['\n-------------------------------------------------------------' repmat('------------', 1, N_samples) ' - \n']);
    fprintf('DIFFERENCE \t%f\t\t\t\t\t\t\t\t\t\t\t', sum(P_G) - (sum(m.P_D(t,:))-wind.P_wf(t)));
    
    difference = sum(R) + wind.P_m(t, :);
    for i = format_samples
        fprintf('%f\t', difference(i));
    end

    fprintf('\n\n\t\t\t\t\t\t\t\tAVERAGE ABSOLUTE DIFFERENCE : \t%f\n\n', meanabs(difference));
    rank_error = 0;
    if isa(m, 'AC_model')
        % print rank Wf
        rank_wf = svd_rank(Wf_opt, 1e2);  % a gap of 100x 
        if rank_wf == 1
            %fprintf('\nRank W_f \t : %i\n', rank_wf);
        else
            rank_error = 1;
            fprintf('\nRank W_f \t : %i\t(!)\n', rank_wf);
        end


        % print rank of each Ws
        for i = format_samples
            rank_ws = svd_rank(Ws_opt(:,:,i), 1e2); % a gap of 100x
            if rank_ws == 1
                %fprintf('Rank W_s %3i : %i\n', i, rank_ws);
            else
                rank_error = 1;
                fprintf('Rank W_s %3i : %i\t(!)\n', i, rank_ws);
            end
        end

        if ~rank_error
            fprintf('All tested W`s are rank 1\n')
        end
    end

   
    fprintf('\nPlotting.');
    h = figure(3);
    set(h, 'name', 'Generators');
    clf
    dock
    fprintf('.');

    plotdata_tmp = [P_G R(:, format_samples)];
    % filter out empty rows
    nonzeros = find(max(plotdata_tmp, [], 2) > 1e-3);
    if length(nonzeros) == 1
        plotdata = [plotdata_tmp(nonzeros, :); nan(1, N+1)];
        labels = {num2str(m.Gens(nonzeros)), ''};
    else
        plotdata = plotdata_tmp(nonzeros, :);
        labels = num2str(m.Gens(nonzeros));
    end

    h = bar(plotdata);
    fprintf('.');

    set(gca, 'xticklabels', labels);
    xlabel('Generator');
    ylabel('Power output');
    grid on;
    title('Generator dispatch');
    legend('P_G', 'R');
    fprintf(repmat('\b', 1, 11));

end

