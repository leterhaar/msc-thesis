%% check the inequality constraint checker 
addpath('../experiment');
N = 100;
init_experiment('model_name', 'case14a', 'wind_N', N);
% define decision variable
x = sdpvar(5*dc.N_G, 1);

% loop over constraint types
N_j = 4*dc.N_G + 2*dc.N_l;
C_ineqs = [];
for j = 1:N_j

    % loop over scenarios
    for i = 1:N
        
        % define constraints
        C_ineqs = [C_ineqs, DC_f_ineq(x, i, dc, wind, 1, j)];
        
    end
end

Obj = DC_f_obj(x, dc, wind, 1);
% optimize (just find feasible point)
optimize([C_ineqs, DC_f_0(x, dc, wind, 1)], Obj, sdpsettings('verbose', 1))


%% check
% loop over constraint types
for j = 1:N_j

    % loop over scenarios
    for i = 1:N
    
        % define constraint and check
        C = DC_f_ineq(x, i, dc, wind, 1, j);
        res_c = check(C);

        % check function
        [~, res_f] = DC_f_check(x, i, dc, wind, 1, j);

        % compare, output if different
        if not(res_c == res_f)
            C
            [res_c res_f abs(res_c - res_f)]
        end
    end
end
            
    
