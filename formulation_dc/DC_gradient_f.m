function grad = DC_gradient_f(x, dc, wind)
% grad = DC_gradient_f(x, dc, wind)
% returns the value of the objective function

    
    
    N_t = size(wind.P_m, 1);
    
    % reshape x if it is a column vector
    if size(x, 2) == 1
        x = reshape(x, 5*dc.N_G, N_t);
    end
    
    grad = zeros_like(x);
    
    % get indices for x
    Rus_idx = dc.N_G+1:2*dc.N_G;
    Rds_idx = 2*dc.N_G+1:3*dc.N_G;

    % loop over time
    for t = 1:N_t
        
        % loop over generators and add linear and quadratic cost gradient
        for k = 1:dc.N_G
            grad(k, t) = 2*dc.c_qu(k) * x(k, t) + dc.c_li(k);
        end
        
    end
    
    % add reserve cost gradient
    grad(Rus_idx, :) = repmat(dc.c_us, 1, N_t);
    grad(Rds_idx, :) = repmat(dc.c_ds, 1, N_t);
end
