function Obj = AC_f_obj(W_f, R_us, R_ds, ac, wind, t, alpha)
% returns the objective function
% 
% Parameters
% ==========
% W_f : forecasted state of network
% R_us : upspinning reserve bounds
% R_ds : downspinning reserve bounds
% ac : instance of AC_model
% wind : instance of wind_model
% t : current time step
% alpha : (optional) if this is given, the epigraph notation will be used
%           for the cost of generator power
%
% Returns
% =======
% Obj = sdpvar (scalar) with objective function

    Obj = 0;
    lambda = 1;

    
    if nargin < 5
        % add power generation costs
        for j = 1:ac.N_G
            k = ac.Gens(j);
            Obj = Obj + ac.c_us(j)*(trace(ac.Y_k(k)*W_f) + ac.P_D(t,k) ...
                                                    -ac.C_w(k)*wind.P_wf(t));
        end
    else
        Obj = sum(alpha);
    end
    
    % add reserve requirements costs
    Obj = Obj + lambda * (ac.c_us' * R_us + ac.c_ds' * R_ds);
end