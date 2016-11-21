function problem = AC_check(x, ac, wind, t)
% tests the solution for feasibility and rank constraints
% returns
% --------
%
% 0 : all good
% 1 : Wf is not PSD
% 2 : Ws is not PSD
% 3 : some W is not rank 1
% 4 : deterministic constraints violated
% 5 : scenario constraints violated

    % extract variables from x
    Wf = x{1};
    Wmus = zero_for_nan(x{2});
    Wmds = zero_for_nan(x{3});

    % check if Wf is psd
    if not(is_psd(Wf))
        problem = 1;
        return
    end

    % check if Wf has rank 1
    if not(svd_rank(Wf, 1e2) == 1)
        problem = 3;
        return
    end

    % create Ws
    N = size(wind.P_w, 2);
    for i = 1:N
        Ws = Wf + Wmus * max(0, -wind.P_m(t, i)) - Wmds * max(0, wind.P_m(t,i));
        if not(svd_rank(Ws, 1e2) == 1)
            problem = 3;
            return
        end
    end

    % test infeasibility deterministic constraints
    x_sdp = {   sdpvar(2*ac.N_b), ...       Wf
        sdpvar(2*ac.N_b), ...       Wmus
        sdpvar(2*ac.N_b), ...       Wmds
        sdpvar(2*ac.N_G, 1)}; ...   Rus and Rds
    C_det = AC_cons_det(x_sdp, ac, wind, t);
    assign_cell(x_sdp, x);
    residuals_det = check(C_det);
    
    if any(residuals_det < -1e-6)
        problem = 4;
        return
    end
    
    % test infeasibility scenario constraints
    if N > 1
        for i = 1:N
            residuals = AC_g(x, ac, wind.slice(i), t);
            if any(residuals > 1e-6)
                problem = 5;
                return
            end
        end
    else
        residuals = AC_g(x, ac, wind, t);
        if any(residuals > 1e-6)
            problem = 5;
            return
        end
    end
        

    problem = 0;
end