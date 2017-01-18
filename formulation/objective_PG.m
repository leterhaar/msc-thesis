function Obj = objective_PG(W_f)
% Obj = objective_PG(W_f)
%
% Evaluates the deterministic part of the objective function
% gets ac model and wind and t from the caller

    check_class({W_f}, {'sdpvar|double'});
    ac = evalin('caller', 'ac');
    wind = evalin('caller', 'wind');
    t = evalin('caller', 't');
    Obj = 0;
    for j = 1:ac.N_G
        k = ac.Gens(j);
        P_Gk = trace(ac.Y_k(k)*W_f) + ac.P_D(t,k)-ac.C_w(k)*wind.P_wf(t);
        Obj = Obj + ac.c_li(j)*P_Gk + ac.c_qu(j) * P_Gk^2;
    end
end