function [C, params] = AC_f_scen_P1(x, ac, wind, t, i, p_des)
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
% i : the scenario index
% p : (optional) the constraint index. If left out, all the constriants
%     corresponding to scenario i will be returned
% 
% Returns
% =======
% C : LMIs with all the constraints
    if nargin < 6
        p_des = 0;
    end

    % extract variables from x
    assert(isa(x, 'cell'), 'x must be cell!');
    n_x = length(x);
    assert(n_x >= 4 && n_x <= 6, 'x must have between 4 and 6 entries');
    W_f = x{1};
    R_us = x{2};
    R_ds = x{3};

    C = [];
    params = zeros(6*ac.N_b+2*ac.N_G+1, 2);
    p = 1;
        
    for k = 1:ac.N_b
        
        % P_inj lower (1')
        if p == p_des || p_des == 0
            C = [C, (ac.P_min(k) ...
                <= trace(ac.Y_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                         + ac.P_D(t, k) - ac.C_w(k)*wind.P_w(t, i)):...
                 sprintf('P_inj lower b%i s%i',k,i)];
            params(p,:) = [i p];
        end
        p = p + 1;
        
        % P_inj upper (1')
        if p == p_des || p_des == 0
            C = [C, (trace(ac.Y_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                         + ac.P_D(t, k) - ac.C_w(k)*wind.P_w(t, i) <= ...
                   ac.P_max(k)):sprintf('P_inj upper b%i s%i',k,i)];
            params(p,:) = [i p];
        end
        p = p + 1;

        % Q_inj lower (2')
        if p == p_des || p_des == 0
            C = [C, (ac.Q_min(k) ...
                 <= trace(ac.Ybar_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                            + ac.Q_D(t, k)):...
                sprintf('Q_inj lower b%i s%i',k,i)];
            params(p,:) = [i p];
        end
        p = p + 1;
        
        % Q_inj upper (2')
        if p == p_des || p_des == 0
            C = [C, (trace(ac.Ybar_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                            + ac.Q_D(t, k) <= ...
                    ac.Q_max(k)):sprintf('Q_inj upper b%i s%i',k,i)];
            params(p,:) = [i p];
        end
        p = p + 1;


        % V_bus lower (3')
        if p == p_des || p_des == 0
            C = [C, (ac.V_min(k)^2 ...
                 <= trace(ac.M_k(k)*(W_f + W_m * wind.P_m(t, i)))):...
                 sprintf('V_bus lower b%i s%i',k,i)];
            params(p,:) = [i p];
        end
        p = p + 1;   
        
        % V_bus upper (3')
        if p == p_des || p_des == 0
            C = [C, (trace(ac.M_k(k)*(W_f + W_m * wind.P_m(t, i))) <= ...
                    ac.V_max(k)^2):sprintf('V_bus upper b%i s%i',k,i)];
            params(p,:) = [i p];
        end
        p = p + 1;   

    end
        

    for g = 1:ac.N_G

        % bus index
        k = ac.Gens(g);

        % Lower bound R between R_us and R_ds
        if p == p_des || p_des == 0
            C = [C, (-R_ds(g) ...
                 <= trace(ac.Y_k(k)*(W_m*wind.P_m(t, i))) ...
                                    - ac.C_w(k)*wind.P_m(t, i)):...
                 sprintf('R lower g%i s%i',g,i)];
            params(p,:) = [i p];
        end
        p = p + 1;
        
        % Lower bound R between R_us and R_ds
        if p == p_des || p_des == 0
            C = [C, (trace(ac.Y_k(k)*(W_m*wind.P_m(t, i))) ...
                                    - ac.C_w(k)*wind.P_m(t, i) <= ...
                    R_us(g)):sprintf('R upper g%i s%i',g,i)];
            params(p,:) = [i p];
        end
        p = p + 1;
        
    end
    
    if p == p_des || p_des == 0
        C = [C, (W_f + W_m * wind.P_m(t, i) >= 0):sprintf('PSD W_s s%i',i)];
        params(p,:) = [i p];
    end
end
