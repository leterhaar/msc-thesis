function [residuals, labels] = DC_g(x, dc, wind, j_des)
%  [residuals, labels] = DC_g(x, dc, wind, j_des)
% Calculates constraint function for one scenario
    
    if nargin < 4
        j_des = 0;
    end
    
    assert(size(wind.P_m, 2) == 1, 'Can only handle one scenario');

    % get indices for x
    PG_idx = 1:dc.N_G;
    Rus_idx = dc.N_G+1:2*dc.N_G;
    Rds_idx = 2*dc.N_G+1:3*dc.N_G;
    dus_idx = 3*dc.N_G+1:4*dc.N_G;
    dds_idx = 4*dc.N_G+1:5*dc.N_G;
       
    % get sizes
    N_j = 4*dc.N_G + 2*dc.N_l;
    N_t = size(wind.P_m, 1);
    
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

        % define reserve power
        R = x(dus_idx, t) * max(0, -wind.P_m(t)) ...
                        - x(dds_idx, t) * max(0, wind.P_m(t));                


        % define scenario power injection vector
        P_injs = dc.C_G * (x(PG_idx, t) + R) ...
                  + dc.C_w * wind.P_w(t) - dc.P_D(t, :)';

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

    
    

