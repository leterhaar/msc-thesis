function Obj = objective(W_f, R_us, R_ds)
%  obj = objective(W_f, R_us, R_ds)

    check_class({W_f, R_us, R_ds}, ...
                {'sdpvar|double', 'sdpvar|double', 'sdpvar|double'});

    Obj = 0;
    
    % get network and wind model and wind time
    ac = evalin('caller', 'ac');
    wind = evalin('caller', 'wind');
    t = evalin('caller', 't');
    

    for j = 1:ac.N_G
            k = ac.Gens(j);
            P_Gk = trace(ac.Y_k(k)*W_f) + ac.P_D(t,k)-ac.C_w(k)*wind.P_wf(t);
            Obj = Obj + ac.c_li(j)*P_Gk + ac.c_qu(j) * P_Gk^2;
    end
     % add reserve requirements costs
    Obj = Obj + ac.c_us' * R_us + ac.c_ds' * R_ds;
end
    
    

   