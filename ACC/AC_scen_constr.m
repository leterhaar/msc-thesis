function C = AC_scen_constr(x, ac, wind, t_wind, scen_id, j_des)
% C = AC_scen_constr(x, ac, wind, t_wind, j_des)
% creates the inequality constraints

    % get the parts between the boxes
    g = AC_g(x, ac, wind, t_wind, j_des);

    if nargin < 6
        j_des = 0;
    end    
    
    C = [];
    j = 1;
        
    for k = 1:ac.N_b
        
        % P_inj lower (1')
        if j == j_des || j_des == 0
            C = [C, (ac.P_min(k) <= g(j)):...
                 sprintf('P_inj lower b%i s%i',k,scen_id)];
        end
        j = j + 1;
        
        % P_inj upper (1')
        if j == j_des || j_des == 0
            C = [C, (g(j) <= ac.P_max(k)):...
                sprintf('P_inj upper b%i s%i',k,scen_id)];
        end
        j = j + 1;

        % Q_inj lower (2')
        if j == j_des || j_des == 0
            C = [C, (ac.Q_min(k) ...
                 <= g(j)):...
                sprintf('Q_inj lower b%i s%i',k,scen_id)];
            params(j,:) = [scen_id j];
        end
        j = j + 1;
        
        % Q_inj upper (2')
        if j == j_des || j_des == 0
            C = [C, (trace(ac.Ybar_k(k)*(W_f + W_m * wind.P_m(t, scen_id))) ...
                            + ac.Q_D(t, k) <= ...
                    ac.Q_max(k)):sprintf('Q_inj upper b%i s%i',k,scen_id)];
            params(j,:) = [scen_id j];
        end
        j = j + 1;


        % V_bus lower (3')
        if j == j_des || j_des == 0
            C = [C, (ac.V_min(k)^2 ...
                 <= trace(ac.M_k(k)*(W_f + W_m * wind.P_m(t, scen_id)))):...
                 sprintf('V_bus lower b%i s%i',k,scen_id)];
            params(j,:) = [scen_id j];
        end
        j = j + 1;   
        
        % V_bus upper (3')
        if j == j_des || j_des == 0
            C = [C, (trace(ac.M_k(k)*(W_f + W_m * wind.P_m(t, scen_id))) <= ...
                    ac.V_max(k)^2):sprintf('V_bus upper b%i s%i',k,scen_id)];
            params(j,:) = [scen_id j];
        end
        j = j + 1;   

    end
        

    for g = 1:ac.N_G

        % bus index
        k = ac.Gens(g);

        % Lower bound R between R_us and R_ds
        if j == j_des || j_des == 0
            C = [C, (-R_ds(g) ...
                 <= trace(ac.Y_k(k)*(W_m*wind.P_m(t, scen_id))) ...
                                    - ac.C_w(k)*wind.P_m(t, scen_id)):...
                 sprintf('R lower g%i s%i',g,scen_id)];
            params(j,:) = [scen_id j];
        end
        j = j + 1;
        
        % Lower bound R between R_us and R_ds
        if j == j_des || j_des == 0
            C = [C, (trace(ac.Y_k(k)*(W_m*wind.P_m(t, scen_id))) ...
                                    - ac.C_w(k)*wind.P_m(t, scen_id) <= ...
                    R_us(g)):sprintf('R upper g%i s%i',g,scen_id)];
            params(j,:) = [scen_id j];
        end
        j = j + 1;
        
    end
    
    if j == j_des || j_des == 0
        C = [C, (W_f + W_m * wind.P_m(t, scen_id) >= 0):sprintf('PSD W_s s%i',scen_id)];
        params(j,:) = [scen_id j];
    end
end


end