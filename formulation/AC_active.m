function params = AC_active(x, ac, wind, t_wind, j_des)
% params = AC_active(x, ac, wind, t_wind, j_des)
% returns the parameters of the active constraints
    
    % tolerance
    tol = 1e-3;
    
    if nargin < 5
        j_des = 0;
    end
    
    params = [];
    
    % only check the psd constraint
    if j_des < 0
        [~, Ws] = AC_g(x, ac, wind, t_wind, 1);
        
        smallest_eig = min(eig(Ws));
        if smallest_eig < tol && smallest_eig > -tol
            params = -1;
        end
    % check for a specific j
    elseif j_des > 0 
        g = AC_g(x, ac, wind, t_wind, j_des);
        if g < tol && g > -tol
            params = j_des;
        end
    % check for all constraints belonging to a scenario
    else
        [g, Ws] = AC_g(x, ac, wind, t_wind);
        all_params = [[1:length(g)]'; -1];
        residuals = [g; min(eig(Ws))];
        params = all_params(residuals < tol & residuals > -tol);
    end
end