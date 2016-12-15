function Obj = DC_f(x, dc, wind)
% Obj = DC_f(x, dc, wind)
% returns the value of the objective function

    N_t = size(wind.P_m, 1);
    
    % reshape x if it is a column vector
    if size(x, 2) == 1
        x = reshape(x, 5*dc.N_G, N_t);
    end
    
    Obj = 0;

    
    % loop over time
    for t = 1:N_t
        % loop over generators and add generator cost at time t
        for k = 1:dc.N_G
            Obj = Obj + (dc.c_qu(k) * (x(k, t))^2) + ...
                                    (dc.c_li(k) * x(k, t));
        end

        % add reserve requirements costs
        Rus = dc.N_G+1:2*dc.N_G;
        Rds = 2*dc.N_G+1:3*dc.N_G;
        Obj = Obj + (dc.c_us' * x(Rus, t) + dc.c_ds' * x(Rds, t));
    end
end
