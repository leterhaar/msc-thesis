function [residuals, labels] = DC_g_delta(x, dc, delta, j_des)
%  [residuals, labels] =  DC_g_delta(x, dc, delta, j_des)
% Calculates residuals for a delta
% delta has the form [P_w1 ... PwN P_mp1 ... P_mpN P_mn1 ... P_mnN]
% where P_mp = max(0, P_m) and P_mn = max(0, -P_m)
    
    if nargin < 4
        j_des = 0;
    end
    
    % get indices for x
    PG_idx = 1:dc.N_G;
    Rus_idx = dc.N_G+1:2*dc.N_G;
    Rds_idx = 2*dc.N_G+1:3*dc.N_G;
    dus_idx = 3*dc.N_G+1:4*dc.N_G;
    dds_idx = 4*dc.N_G+1:5*dc.N_G;

    % get sizes
    N_j = 4*dc.N_G + 2*dc.N_l;
    N_t = size(delta, 2)/3;

    % reshape x if it is a column vector
    if size(x, 2) == 1
        x = reshape(x, 5*dc.N_G, N_t);
    end

    % preallocate residual vector
    if isa(x, 'sdpvar')
        residuals = sdpvar(N_t*N_j, 1);
    else
        residuals = nan(N_t*N_j, 1);
        x = zero_for_nan(x); % remove nans
    end
    labels = cell(N_t*N_j, 1);
    
    % loop over time steps
    j = 1;
    for t = 1:N_t

        % define reserve power (negative mismatch starts from 2N_tt, positive
        % from N_t
                            
        R = x(dus_idx, t) * delta(2*N_t+t) ...
                        - x(dds_idx, t) * delta(N_t+t);                


        % define scenario power injection vector
        P_injs = dc.C_G * (x(PG_idx, t) + R) ...
                  + dc.C_w * delta(t) - dc.P_D(t, :)';

        for k = 1:dc.N_G   

            % generator lower limit
            if j_des == j || j_des == 0
                residuals(j) = - dc.P_Gmin(k) + x(PG_idx(k), t) + R(k);
                labels{j} = sprintf('Gen lower g%2i t%2i', k, t);
            end
            j = j+1;

            % generator upper limit
            if j_des == j || j_des == 0
                residuals(j) = - x(PG_idx(k), t) - R(k) + dc.P_Gmax(k);
                labels{j} = sprintf('Gen upper g%2i t%2i', k, t);
            end
            j = j+1;    

            % R lower bound
            if j_des == j || j_des == 0
                residuals(j) = x(Rds_idx(k), t) + R(k);
                labels{j} = sprintf('R lower g%2i t%2i', k, t);
            end
            j = j+1;   

            % R upper bound
            if j_des == j || j_des == 0
                residuals(j) = - R(k) + x(Rus_idx(k), t);
                labels{j} = sprintf('R upper g%2i t%2i', k, t);
            end
            j = j+1;   

        end

        P_fs = dc.B_f * [dc.B_bustildeinv * P_injs(1:end-1); 0];

        for k = 1:dc.N_l

            % line flow lower limit
            if j_des == j || j_des == 0
                residuals(j) = dc.P_fmax(k) + P_fs(k);
                labels{j} = sprintf('Pf lower l%2i t%2i', k, t);
            end
            j = j+1;   

            % line flow lower limit
            if j_des == j || j_des == 0
                residuals(j) = - P_fs(k) + dc.P_fmax(k);
                labels{j} = sprintf('Pf upper l%2i t%2i', k, t);
            end
            j = j+1;   

        end
    end
    
    if j_des > 0
        residuals = residuals(j_des);
    end
end

    
    

