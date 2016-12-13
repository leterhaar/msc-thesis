function C = DC_cons_scen_delta(x, dc, delta, j_des)
% C = DC_cons_scen_delta(x, dc, delta, j_des)
% returns an LMI with all the constraints, not aggregated
% if j_des < 0, returns only the psd constraint

    if nargin < 4
        j_des = 0;
    end
    
    C = [];
    
    % retrieve 1
    if j_des > 0
        [g, labels] = DC_g_delta(x, dc, delta, j_des);
        C = [(g >= 0):labels{j_des}];
    
    % retrieve all
    else
        [g, labels] = DC_g_delta(x, dc, delta);
        
        for j = 1:length(g)
            C = [C, (g(j) >= 0):labels{j}];
        end
    end
    
end