function [C, params] = AC_f_ineq(x, i, ac, wind, t, j_des)
% [C, params] = DC_f_ineq(x, i, dc, wind, t_wind)
% Function to create a constraint set based on a scenario

    if nargin < 6
        j_des = 0;
    end

    % extract variables from x
    W_f = x{1};
    W_m = x{2};
    R = x{3};
    R_us = R(1:ac.N_G);
    R_ds = R(ac.N_G+1:2*ac.N_G);
    
    
    C = [];
    params = zeros(6*ac.N_b+2*ac.N_G+1, 2);
    j = 1;
        
    for k = 1:ac.N_b
        
        % P_inj lower (1')
        if j == j_des || j_des == 0
            C = [C, (ac.P_min(k) ...
                <= trace(ac.Y_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                         + ac.P_D(t, k) - ac.C_w(k)*wind.P_w(t, i)):...
                 sprintf('P_inj lower b%i s%i',k,i)];
            params(j,:) = [i j];
        end
        j = j + 1;
        
        % P_inj upper (1')
        if j == j_des || j_des == 0
            C = [C, (trace(ac.Y_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                         + ac.P_D(t, k) - ac.C_w(k)*wind.P_w(t, i) <= ...
                   ac.P_max(k)):sprintf('P_inj upper b%i s%i',k,i)];
            params(j,:) = [i j];
        end
        j = j + 1;

        % Q_inj lower (2')
        if j == j_des || j_des == 0
            C = [C, (ac.Q_min(k) ...
                 <= trace(ac.Ybar_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                            + ac.Q_D(t, k)):...
                sprintf('Q_inj lower b%i s%i',k,i)];
            params(j,:) = [i j];
        end
        j = j + 1;
        
        % Q_inj upper (2')
        if j == j_des || j_des == 0
            C = [C, (trace(ac.Ybar_k(k)*(W_f + W_m * wind.P_m(t, i))) ...
                            + ac.Q_D(t, k) <= ...
                    ac.Q_max(k)):sprintf('Q_inj upper b%i s%i',k,i)];
            params(j,:) = [i j];
        end
        j = j + 1;


        % V_bus lower (3')
        if j == j_des || j_des == 0
            C = [C, (ac.V_min(k)^2 ...
                 <= trace(ac.M_k(k)*(W_f + W_m * wind.P_m(t, i)))):...
                 sprintf('V_bus lower b%i s%i',k,i)];
            params(j,:) = [i j];
        end
        j = j + 1;   
        
        % V_bus upper (3')
        if j == j_des || j_des == 0
            C = [C, (trace(ac.M_k(k)*(W_f + W_m * wind.P_m(t, i))) <= ...
                    ac.V_max(k)^2):sprintf('V_bus upper b%i s%i',k,i)];
            params(j,:) = [i j];
        end
        j = j + 1;   

    end
        

    for g = 1:ac.N_G

        % bus index
        k = ac.Gens(g);

        % Lower bound R between R_us and R_ds
        if j == j_des || j_des == 0
            C = [C, (-R_ds(g) ...
                 <= trace(ac.Y_k(k)*(W_m*wind.P_m(t, i))) ...
                                    - ac.C_w(k)*wind.P_m(t, i)):...
                 sprintf('R lower g%i s%i',g,i)];
            params(j,:) = [i j];
        end
        j = j + 1;
        
        % Lower bound R between R_us and R_ds
        if j == j_des || j_des == 0
            C = [C, (trace(ac.Y_k(k)*(W_m*wind.P_m(t, i))) ...
                                    - ac.C_w(k)*wind.P_m(t, i) <= ...
                    R_us(g)):sprintf('R upper g%i s%i',g,i)];
            params(j,:) = [i j];
        end
        j = j + 1;
        
    end
    
    if j == j_des || j_des == 0
        C = [C, (W_f + W_m * wind.P_m(t, i) >= 0):sprintf('PSD W_s s%i',i)];
        params(j,:) = [i j];
    end
end
