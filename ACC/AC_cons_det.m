function C = AC_cons_det(x, ac, wind, t)
% returns the deterministic constraints

    % extract variables from x
    W_f = x{1};
    W_mus = x{2};
    W_mds = x{3};
    R = x{4};

    C = [];
    k_ref = ac.N_b + ac.refbus;

    for k = 1:ac.N_b
        % P_inj (1)
        C = [C, (ac.P_min(k) ...
            <= trace(ac.Y_k(k)*W_f) + ac.P_D(t, k) - ac.C_w(k)*wind.P_wf(t) <= ...
               ac.P_max(k)):...
               sprintf('P_inj b%i det', k)];

        % Q_inj (2)
        C = [C, (ac.Q_min(k) ...
            <= trace(ac.Ybar_k(k)*W_f) + ac.Q_D(t, k) <=...
               ac.Q_max(k)):...
               sprintf('Q_inj b%i det', k)];

        % V_bus (3)
        C = [C, (ac.V_min(k)^2 ...
            <= trace(ac.M_k(k)*W_f) <= ...
               ac.V_max(k)^2):...
               sprintf('V_bus b%i det', k)];
    end
    % PSD constraint on W_f 
    C = [C, (W_f >= 0):'PSD on Wf'];

    % refbus constraint
    C = [C, (W_f(k_ref, k_ref) == 0):'Refbus Wf'];
    C = [C, (W_mus(k_ref, k_ref) == 0):'Refbus Wmus'];
    C = [C, (W_mds(k_ref, k_ref) == 0):'Refbus Wmds'];

    Ysum = zeros_like(ac.Y_k(1));
    for k = ac.Gens'
        Ysum = Ysum + ac.Y_k(k);
    end

    % sum Wm = 1
    C = [C, (trace(Ysum * W_mus) == 1):'Balancing Wmus'];   
    C = [C, (trace(Ysum * W_mds) == 1):'Balancing Wmds'];

    % Nonnegativity constraints on reserve bounds
    C = [C, (R >= 0):'NonNeg R'];
end