function C = AC_f_det(x, ac, wind, t)
% returns the deterministic constraints
%
% Parameters
% ==========
% x : a cell with the following entries:
%       - W_f : forecasted state of network
%       - R_us : upspinning reserve bound
%       - R_ds : downspinning reserve bound
%       Using P1: 
%       - W_s : a [2Nb x 2Nb x N] matrix with the scenario network states
%       - d_us : the [NG x 1] upspinning distribution vector
%       - d_ds : the [NG x 1] downspinning distribution vector
%       Using P2:
%       - W_m, a [2Nb x 2Nb] matrix with the mismatch network 
%         state
%       Using P3:
%       - W_mus, a [2Nb x 2Nb] matrix with the upspinning 
%         mismatch network state
%       - W_mds, a [2Nb x 2Nb] matrix with the downspinning
%         mismatch network state
% ac : instance of AC_model
% wind : instance of wind_model
% t : current time step
% 
% Returns
% =======
% C : LMIs with all the constraints

    % extract variables from x
    assert(isa(x, 'cell'), 'x must be cell!');
    n_x = length(x);
    assert(n_x >= 4 && n_x <= 7, 'x must have between 4 and 6 entries');
    W_f = x{1};
    R_us = x{2};
    R_ds = x{3};
    
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
    C = [C, (W_f >= 0):'PSD on W_f'];
    
    % Nonnegativity constraints on reserve bounds
    C = [C, (R_us >= 0):'NN Rus', (R_ds >= 0):'NN Rds'];

    % refbus constraint
    C = [C, (W_f(k_ref, k_ref) == 0):'Refbus Wf'];
    
    % FORMULATION SPECIFIC CONSTRAINTS
    
    % P1
    if n_x == 6
        % x{4} is W_s, not used in deterministic constraints
        d_us = x{5};
        d_ds = x{6};
        
        C = [C, (sum(d_us) == 1):'Sum us'];
        C = [C, (sum(d_ds) == 1):'Sum ds'];
    
    % P2
    elseif n_x == 4
        W_m = x{4};
        
        Ysum = zeros_like(ac.Y_k(1));
        for k = ac.Gens'
            Ysum = Ysum + ac.Y_k(k);
        end

        C = [C, (W_m(k_ref, k_ref) == 0):'Refbus Wm'];
        % sum Wm = 1
        C = [C, (trace(Ysum * W_m) == -1):'Balancing constraints'];
        
    % P3
    else
        W_mus = x{4};
        W_mds = x{5};
        
        Ysum = zeros_like(ac.Y_k(1));
        for k = ac.Gens'
            Ysum = Ysum + ac.Y_k(k);
        end

        C = [C, (W_mus(k_ref, k_ref) == 0):'Refbus Wmus'];
        C = [C, (W_mds(k_ref, k_ref) == 0):'Refbus Wmds'];
        
        C = [C, (trace(Ysum * W_mus) == 1):'Balancing constraints us'];   
        C = [C, (trace(Ysum * W_mds) == 1):'Balancing constraints ds'];
    end
end