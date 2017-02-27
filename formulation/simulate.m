function [violations, amount] = simulate(ac, Wf, d_us, d_ds, wind, verbose)
% runs a MC simulation using MATPOWERs runpf
% for one wind trajectory
%
% ARGUMENTS
% =========
% ac : ac model
% Wf : cell with T state matrices for forecasted state
% d_us : cell with T upspinning reserve distribution vectors
% d_ds : cell with T downspinning reserve distribution vectors
% wind : wind model
% verbose : flag to suppress progress, default true

    if nargin < 6
        verbose = 1;
    end
    
    % prepare
    new_mpc = loadcase(ac.model_name);
    [T, N_s] = size(wind.P_w);

    violations = nan(ac.N_l, T);
    amount = nan(ac.N_l, T, N_s+1);
    opts = mpoption('verbose', 0, 'out.all', 0);
    
    if verbose
        prog = progress('Simulating power flows', T*(N_s+1));
    end
    % forecast
    % loop over time steps
    for t = 1:T
        
        % reset mpc
        mpc = new_mpc;
        
        % set loads (PQ buses)
        mpc.bus(:, 3) = (ac.P_D(t,:)' - ac.C_w*wind.P_wf(t))*mpc.baseMVA;
        mpc.bus(:, 4) = ac.Q_D(t,:)*mpc.baseMVA;

        
        % set generators (PV buses)
        [Pg, Qg, Vmag] = extract_dispatch(Wf{t},t);
        for j = 1:ac.N_G
            k = ac.Gens(j);
            mpc.gen(j, 2) = Pg(j)*mpc.baseMVA;
            mpc.gen(j, 3) = Qg(j)*mpc.baseMVA;
            mpc.gen(j, 6) = Vmag(k);
        end
        
        results = runpf(mpc, opts);
        
        % find power flows over lines
        P_from = results.branch(:,14);
        P_to = results.branch(:,16);
        Q_from = results.branch(:,15);
        Q_to = results.branch(:,17);
        
        % define maximum (NB ask Vahab whether this is ok)
        P_f = max(abs([P_from P_to]), [], 2);
        Q_f = max(abs([Q_from Q_to]), [], 2);
        
        % calculate app line flows and compare with limits
        S_f = sqrt(P_f.^2 + Q_f.^2);
        
        % check 
        violations(:, t) = S_f > (ac.S_max*mpc.baseMVA);
        amount(:, t, 1) = S_f / ac.S_max*mpc.baseMVA * 100;
        
        if verbose
            prog.ping();
        end
        
        
        % loop over scenarios and do the same
        for i = 1:N_s
            
            % set wind power to different realization
            mpc.bus(:, 3) = (ac.P_D(t,:)' - ac.C_w*wind.P_w(t,i))*mpc.baseMVA;
            
            % define reserve power
            R = -min(wind.P_m(t,i), 0)*d_us{t} - max(wind.P_m(t,i), 0)*d_ds{t};
            mpc.gen(:, 2) = mpc.gen(:, 2) + R;
            
            results = runpf(mpc, opts);

            % find power flows over lines
            P_from = results.branch(:,14);
            P_to = results.branch(:,16);
            Q_from = results.branch(:,15);
            Q_to = results.branch(:,17);

            % define maximum (NB ask Vahab whether this is ok)
            P_f = max(abs([P_from P_to]), [], 2);
            Q_f = max(abs([Q_from Q_to]), [], 2);

            % calculate app line flows and compare with limits
            S_f = sqrt(P_f.^2 + Q_f.^2);

            % check 
            violations(:, t) = violations(:,t) + (S_f > (ac.S_max*mpc.baseMVA));
            amount(:, t, i+1) = S_f / ac.S_max*mpc.baseMVA * 100;
            
            if verbose
                prog.ping();
            end
        end
    end
    
end