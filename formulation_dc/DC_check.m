function flag = DC_check(x, dc, wind)
% flag = DC_check(x, dc, wind)
% checks the solution against all constraints
% returns
% 0 feasible
% 1 infeasible due to deterministic constraints
% 2 infeasible due to some scenario constraint

    tol = 1e-6;
    [N_t, N] = size(wind.P_m);
    
    % reshape x if it is a column vector
    if size(x, 2) == 1
        x = reshape(x, 5*dc.N_G, N_t);
    end
    
    % check deterministic constraints
    x_sdp = sdpvar(5*dc.N_G, N_t, 'full');
    C_det = DC_cons_det(x_sdp, dc, wind);
    assign(x_sdp, x);
    residuals = check(C_det);
    if any(residuals < -tol)
        flag = 1;
        return
    end
    
    % loop over scenarios and check
    for i = 1:N
        residuals = DC_g(x, dc, wind.slice(i));
        if any(residuals < -tol)
            flag = 2;
            return
        end
    end
    
    % all good
    flag = 0;
end
