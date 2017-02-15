% "Recovery.m" recovers voltage phasors from the matrix solution and
% derives the generating power.

% SDP Solver of Optimal Power Flow: beta version
% Ramtin Madani(madani@ee.columbia.edu)
% Morteza Ashraphijuo(morteza.ashraphijuo@gmail.com)
% Javad Lavaei(lavaei@ee.columbia.edu)
% Columbia University
% Last Modified: October 07, 2014

function [V_rec, Sg_rec, Sb_rec, sf_rec, st_rec, cost_rec] = Recovery(V2,W,Pd,Qd, ...
    ng,incidentG,Qmax,Qmin,Pmax,Pmin,correction,fbusN,tbusN,activeL,busType,Yf,Yt,Ybus,cc, ...
    solver_rec,c2,c1,c0)

    nc = length(cc);
    nb = size(V2,1);
    nl = length(fbusN);

    Vm_rec = sqrt(V2);
    phaseL = zeros(nl,nc);
    for rr = 1 : nc
        in1 = nb*(rr-1);
        in2 = setdiff(activeL, cc{rr});
        phaseL(in2,rr) = angle( diag( W(in1 + fbusN(in2), in1 + tbusN(in2)) ) );
    end
    source = find(busType == 3);

    cvx_begin

    if strcmp(solver_rec,'sdpt3')
        cvx_solver sdpt3
    elseif strcmp(solver_rec,'sedumi')
        cvx_solver sedumi
    elseif strcmp(solver_rec,'mosek')
        cvx_solver mosek
    end

    cvx_precision best

    variable Va_rec(nb,nc)
    expression phaseL1(nl,nc)

    for rr = 1 : nc
        in2 = setdiff(activeL, cc{rr});
        phaseL1(in2,rr) = Va_rec(fbusN(in2),rr) - Va_rec(tbusN(in2),rr);
    end

    cost_lp = sum(sum(abs(phaseL1 - phaseL)));

    minimize( cost_lp );

    subject to

    Va_rec(source,:) == 0;

    cvx_end

    V_rec = Vm_rec .* exp(Va_rec*1i);

    Sb_rec = zeros(nb,nc);
    sf_rec = zeros(nl,nc);
    st_rec = zeros(nl,nc);

    for rr = 1 : nc
        ex = setdiff(1 : nl, cc{rr});    

        Sb_rec(:,rr) = Pd + 1i*Qd + V_rec(:,rr) .* conj(Ybus{rr}*V_rec(:,rr));

        sf_rec(ex,rr) = V_rec(fbusN(ex),rr) .* conj(Yf(ex,:) * V_rec(:,rr));
        st_rec(ex,rr) = V_rec(tbusN(ex),rr) .* conj(Yt(ex,:) * V_rec(:,rr));
    end

    gb = find((sum(incidentG) ~= 0));

    cvx_begin

    if strcmp(solver_rec,'sdpt3')
        cvx_solver sdpt3
    elseif strcmp(solver_rec,'sedumi')
        cvx_solver sedumi
    elseif strcmp(solver_rec,'mosek')
        cvx_solver mosek
    end

    cvx_precision best

    variable Sg_rec(ng,nc) complex
    
    cost_rec = sum(c2.*real(Sg_rec(:,1)).^2 + c1.*real(Sg_rec(:,1)) + c0);

    minimize( cost_rec + (5*10^8) * ...
        max([0;max(max(real(Sg_rec) - Pmax * ones(1,nc))); ...
               max(max(Pmin * ones(1,nc) - real(Sg_rec))); ...
               max(max(imag(Sg_rec) - Qmax * ones(1,nc))); ...
               max(max(Qmin * ones(1,nc) - imag(Sg_rec))); ...
               max(max(abs(Sb_rec(gb,:) - incidentG(:,gb).' * Sg_rec)))]));
           
%     minimize( cost_rec );

    subject to

%     abs(Sb_rec(gb) - incidentG(:,gb).' * Sg_rec) <= 0;

%     real(Sg_rec) <= Pmax * ones(1,nc);
%     real(Sg_rec) >= Pmin * ones(1,nc);
%     imag(Sg_rec) <= Qmax * ones(1,nc);
%     imag(Sg_rec) >= Qmin * ones(1,nc);

    for rr = 2 : nc
        abs(real(Sg_rec(:,rr) - Sg_rec(:,1))) <= correction;
    end    
    
    cvx_end
    
end