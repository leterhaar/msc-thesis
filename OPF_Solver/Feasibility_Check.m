% "Feasibility_Check.m" checks the feasiblity of the derived solution from SDP.

% SDP Solver of Optimal Power Flow: beta version
% Ramtin Madani(madani@ee.columbia.edu)
% Morteza Ashraphijuo(morteza.ashraphijuo@gmail.com)
% Javad Lavaei(lavaei@ee.columbia.edu)
% Columbia University
% Last Modified: October 07, 2014

function [feas_flag, vio] = Feasibility_Check(V_rec,Sg_rec,Sb_rec,sf_rec,st_rec, ...
    incidentG,Vmax,Vmin,activeG,Qmax,Qmin,Pmax,Pmin,activeL,SlmMax,nc,tol_feas,correction)
    
    cpv = zeros(length(activeG),nc-1);
    for rr = 2 : nc
        cpv(:,rr-1) = abs(real(Sg_rec(:,rr) - Sg_rec(:,1)));
    end

    vio.P_balance = sum(abs(real(Sb_rec - incidentG.' * Sg_rec)) >= tol_feas);
    vio.Q_balance = sum(abs(imag(Sb_rec - incidentG.' * Sg_rec)) >= tol_feas);
    vio.gen_Pmax = sum(real(Sg_rec) >= Pmax * ones(1,nc) + tol_feas);
    vio.gen_Pmin = sum(real(Sg_rec) <= Pmin * ones(1,nc) - tol_feas);
    vio.gen_Qmax = sum(imag(Sg_rec) >= Qmax * ones(1,nc) + tol_feas);
    vio.gen_Qmin = sum(imag(Sg_rec) <= Qmin * ones(1,nc) - tol_feas);
    vio.Vmax = sum(abs(V_rec) >= Vmax * ones(1,nc) + tol_feas);
    vio.Vmin = sum(abs(V_rec) <= Vmin * ones(1,nc) - tol_feas);
    vio.SlmMax = sum(abs(sf_rec(activeL,:)) >= SlmMax(activeL) * ones(1,nc) + tol_feas);
    vio.SmlMax = sum(abs(st_rec(activeL,:)) >= SlmMax(activeL) * ones(1,nc) + tol_feas);
    vio.cont = sum(cpv >= correction + tol_feas);
    vio.total = sum(vio.P_balance + vio.Q_balance + vio.gen_Pmax + vio.gen_Pmin + ...
        vio.gen_Qmax + vio.gen_Qmin + vio.Vmax + vio.Vmin + vio.SlmMax + vio.SmlMax) + sum(vio.cont);

    %% Output results
    disp('Active power balance violations:')
    disp(vio.P_balance)
    disp('Reactive power balance violations:')
    disp(vio.Q_balance)
   
    disp('Generator capacity violations Pmax:')
    disp(vio.gen_Pmax)
    disp('Generator capacity violations Pmin:')
    disp(vio.gen_Pmin)
    disp('Generator capacity violations Qmax:')
    disp(vio.gen_Qmax)
    disp('Generator capacity violations Qmin:')
    disp(vio.gen_Qmin)

    disp('Max voltage violations:')
    disp(vio.Vmax)
    disp('Min voltage violations:')
    disp(vio.Vmin)

    disp('Line rating violations:')
    disp(vio.SlmMax)
    disp(vio.SmlMax)

    disp('Controllable parameters violations:')
    disp(vio.cont)
    
    if vio.total == 0
        feas_flag = 1;
    else
        feas_flag = 0;
    end

end
