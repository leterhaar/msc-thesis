function [P_G, Q_G, Vmag] = extract_dispatch(W, t)
% [PGopt, QGopt, Vopt] = extract_dispatch(W, t)
% uses 'wind' and 'ac' from caller to  extract generator dispatch from W at
% time t

    % extract network information from caller
    ac = evalin('caller', 'ac');
    wind = evalin('caller', 'wind');
    
    % find generator dispatch
    P_G = nan(ac.N_G,1);
    Q_G = nan(ac.N_G,1);
    
    for j = 1:ac.N_G
        k = ac.Gens(j);
        P_G(j) = trace(ac.Y_(k)*W) - ac.C_w(k)*wind.P_wf(t) + ac.P_D(t, k);
        Q_G(j) = trace(ac.Y_(k)*W) + ac.Q_D(t, k);
    end
    
    diagonal = diag(W);
    Vmag = sqrt(diagonal(1:ac.N_b) + diagonal(ac.N_b+1:end));
    
end