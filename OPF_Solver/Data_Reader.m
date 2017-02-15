% "Data_Reader.m" reads matpower matpower data file format.

% SDP Solver of Optimal Power Flow: beta version
% Ramtin Madani(madani@ee.columbia.edu)
% Morteza Ashraphijuo(morteza.ashraphijuo@gmail.com)
% Javad Lavaei(lavaei@ee.columbia.edu)
% Columbia University
% Last Modified: October 07, 2014

function [nb,busL,busN,busType,Pd,Qd,Vmax,Vmin,statusG,activeG,ng,genL,genN,incidentG,Qmax,Qmin,...
    Pmax,Pmin,activeL,statusL,nl,fbusL,tbusL,SlmMax,fbusN,tbusN,incidentF,incidentT,...
    Yf,Yt,YfP,YtP,Ybus,edges,c2,c1,c0] = Data_Reader(mpc,cc)
       
    nc = length(cc);

    %Importing bus parameres
    nb = size(mpc.bus,1);
    busL = mpc.bus(:,1);
    busN = sparse(busL, ones(nb,1), (1:nb) );
    busType = mpc.bus(:,2);
    Pd = mpc.bus(:,3)/mpc.baseMVA;
    Qd = mpc.bus(:,4)/mpc.baseMVA;

    Gs = mpc.bus(:,5)/mpc.baseMVA;
    Bs = mpc.bus(:,6)/mpc.baseMVA;
    Vmax = mpc.bus(:,12);
    Vmin = mpc.bus(:,13);

    Ysh = Gs + Bs*1i;

    %Importing generator parameres
    statusG = mpc.gen(:,8);
    activeG = find(statusG);

    ng = size(mpc.gen,1);

    genL = mpc.gen(:,1);
    genN = busN(genL);
    incidentG = diag(statusG)*sparse(1:ng, genN, 1 , ng, nb);

    Qmax = (mpc.gen(:,4).*statusG)/mpc.baseMVA;
    Qmin = (mpc.gen(:,5).*statusG)/mpc.baseMVA;
    Pmax = (mpc.gen(:,9).*statusG)/mpc.baseMVA;
    Pmin = (mpc.gen(:,10).*statusG)/mpc.baseMVA;

    %Importing branch parameres
    statusL = mpc.branch(:,11);
    activeL = find(statusL);

    nl = size(mpc.branch,1);

    fbusL = mpc.branch(:,1);
    tbusL = mpc.branch(:,2);
    rBranch = mpc.branch(:,3);
    xBranch = mpc.branch(:,4);
    bBranch = mpc.branch(:,5);
    SlmMax = mpc.branch(:,6)/mpc.baseMVA;

    ratio = mpc.branch(:,9);
    phase = mpc.branch(:,10);

    ratio(ratio == 0) = 1;

    ys = 1 ./ (rBranch + xBranch*1i);% + 2*10^(-7)*ones(nl,1);
    sh = exp(phase*(pi/180)*1i);

    fbusN = busN(fbusL);
    tbusN = busN(tbusL);

    incidentF = diag(statusL) * sparse(1 : nl, fbusN, 1, nl, nb);
    incidentT = diag(statusL) * sparse(1 : nl, tbusN, 1, nl, nb);

    Yff = (ys + (bBranch/2)*1i) ./ (ratio.^2);
    Yft = -(ys ./ ratio) .* sh;
    Ytf = -(ys ./ ratio) ./ sh;
    Ytt = (ys + (bBranch/2) * 1i);

    Yf = sparse(diag(Yff) * incidentF + diag(Yft) * incidentT);
    Yt = sparse(diag(Ytf) * incidentF + diag(Ytt) * incidentT);

    YffP = ys ./ (ratio .^ 2);             YftP = -ys ./ (ratio .* sh);
    YtfP = -ys ./ (ratio ./ sh);           YttP = ys;
    
    YfP = sparse(diag(YffP) * incidentF + diag(YftP) * incidentT);
    YtP = sparse(diag(YtfP) * incidentF + diag(YttP) * incidentT);

    Ybus = cell(nc,1);
    edges = zeros(nl,nc);
    for rr = 1 : nc
        ex = setdiff(1 : nl, cc{rr});
        ex2 = setdiff(activeL, cc{rr});
        edges(ex2,rr) = 1;
        Ybus{rr} = incidentF(ex, :)' * Yf(ex, :) + incidentT(ex, :)' * Yt(ex, :) + diag(Ysh);
    end

    %Importing cost parameres
    c2 = (mpc.gencost(:,5).*statusG)*mpc.baseMVA^2;
    c1 = (mpc.gencost(:,6).*statusG)*mpc.baseMVA;
    c0 = (mpc.gencost(:,7).*statusG);
    
end
