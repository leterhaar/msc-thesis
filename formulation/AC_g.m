function [g, W_s, labels] = AC_g(x, ac, wind, t, j_des)
% [g, W_s] = AC_g_scen(x, ac, P_wf, P_w, t, j_des)
% formulation of the constraint functions for one scenario g(x,delta) <= 0

    W_f = x{1};
    W_mus = x{2};
    W_mds = x{3};
    R = x{4};
    R_us = R(1:ac.N_G);
    R_ds = R(ac.N_G+1:2*ac.N_G);

    if nargin < 5
        j_des = 0;
    end
    
    assert(size(wind.P_m, 2) == 1, 'Can only handle 1 scenario');

    % preallocate constraints
    N_cons = 6*ac.N_b + 2*ac.N_G;
    if isa(W_f, 'sdpvar')
        g = sdpvar(N_cons, 1);
    else
        g = nan(N_cons, 1);
        % set nans to zero if needed
        W_mus = zero_for_nan(W_mus);
        W_mds = zero_for_nan(W_mds);
    end
    labels = cell(N_cons, 1);
    
    % define scenario state
    W_s = W_f + W_mus * max(0, -wind.P_m(t)) ...
              - W_mds * max(0, wind.P_m(t));
    
    
    % loop over buses
    j = 1;
    for k = 1:ac.N_b
        
        % P_inj upper (1)
        if j_des == j || j_des == 0
            g(j) = ... P_G
                   trace(ac.Y_k(k)*W_s) + ac.P_D(t, k) - ac.C_w(k)*wind.P_w(t) ...
                   ... <= P_Gmax
                   - ac.P_max(k);
            labels{j} = sprintf('P_inj upper b%2i', k);
        end
        j = j+1;
        
        % P_inj lower (1)
        if j_des == j || j_des == 0
            g(j) = ... P_G
                   -(trace(ac.Y_k(k)*W_s) + ac.P_D(t, k) - ac.C_w(k)*wind.P_w(t)) ...
                   ... >= P_Gmin
                   + ac.P_min(k);
            labels{j} = sprintf('P_inj lower b%2i', k);
        end
        j = j+1;

        % Q_inj upper (2) 
        if j_des == j || j_des == 0
            g(j) = ... Q_G
                   trace(ac.Ybar_k(k)*W_s) + ac.Q_D(t, k) ...
                   ... <= Q_Gmax
                   - ac.Q_max(k);
            labels{j} = sprintf('Q_inj upper b%2i', k);
        end
        j = j+1;
        
        % Q_inj lower (2) 
        if j_des == j || j_des == 0
            g(j) = ... Q_G
                   -(trace(ac.Ybar_k(k)*W_s) + ac.Q_D(t, k)) ...
                   ... >= Q_Gmin
                   + ac.Q_min(k);
            labels{j} = sprintf('Q_inj lower b%2i', k);
        end
        j = j+1;

        % V_bus upper (3)
        if j_des == j || j_des == 0
            g(j) = ... V_bus
                   trace(ac.M_k(k)*W_s) ...
                   ... <= (V_max)^2
                   - (ac.V_max(k))^2;
            labels{j} = sprintf('V_bus upper b%2i', k);
        end
        j = j+1;
        
        % V_bus lower (3)
        if j_des == j || j_des == 0
            g(j) = ... V_bus
                   -(trace(ac.M_k(k)*W_s)) ...
                   ... >= (V_min)^2
                   + (ac.V_min(k))^2;
            labels{j} = sprintf('V_bus lower b%2i', k);
        end
        j = j+1;
    end

    % build R cosntraints
    for i = 1:ac.N_G

        % bus index
        k = ac.Gens(i);
      
        % R upper
        if j_des == j || j_des == 0
            g(j) = ... R
                   trace(ac.Y_k(k)*(W_s-W_f)) - ac.C_w(k)*wind.P_m(t) ...
                   ... <= R_us
                   - R_us(i);
            labels{j} = sprintf('R upper g%2i', k);
        end
        j = j+1;
        
        % R lower
        if j_des == j || j_des == 0
            g(j) = ... R
                   - (trace(ac.Y_k(k)*(W_s-W_f)) - ac.C_w(k)*wind.P_m(t)) ...
                   ... >= -R_ds
                   - R_ds(i);
            labels{j} = sprintf('R lower g%2i', k);
        end
        j = j+1;

    end
    if j_des > 0 
        g = g(j_des);
    end
end

