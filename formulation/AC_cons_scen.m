function C = AC_cons_scen(x, ac, wind, t_wind, j_des)
% C = AC_cons_scen(x, ac, wind, t_wind, j_des)
% returns an LMI with all the constraints
% if j_des < 0, returns only the psd constraint

    if nargin < 5
        j_des = 0;
    end
    
    C = [];
    
    % retrieve psd constraint
    if j_des < 0
        [~, Ws] = AC_g(x, ac, wind, t_wind, 1);
        C = [(Ws >= 0):'Ws >= 0'];
        
    % retrieve specific constraint
    elseif j_des > 0
        [g, ~, labels] = AC_g(x, ac, wind, t_wind, j_des);
        C = [(g <= 0):labels{j_des}];
    
    % retrieve all
    else
        [g, Ws, labels] = AC_g(x, ac, wind, t_wind);
        
        for j = 1:length(g)
            C = [C, (g(j) <= 0):labels{j}];
        end
        
        % add psd
        C = [C, [Ws >= 0]:'Ws >= 0'];
    end
    
end