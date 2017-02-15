% "OPF_Solver.m" solves the opf problem given the input parameters.

% SDP Solver of Optimal Power Flow: beta version
% Ramtin Madani(madani@ee.columbia.edu)
% Morteza Ashraphijuo(morteza.ashraphijuo@gmail.com)
% Javad Lavaei(lavaei@ee.columbia.edu)
% Columbia University
% Last Modified: October 07, 2014

function [results] = OPF_Solver(data,varargin)

    warning off
    
    mpc = eval(data);

    if nargin == 2
	    opt_set = varargin{1};
        if isfield(opt_set,'epB')
            ep(1) = opt_set.epB;
        else
            ep(1) = 0;
        end
        if isfield(opt_set,'epL')
            ep(2) = opt_set.epL;
        else
            ep(2) = 0;
        end
        if isfield(opt_set,'cc')
            cc = opt_set.cc;
        else
            cc{1} = double.empty(1,0);
        end
        ndxL_prob = zeros(length(mpc.branch(:,1)), size(cc,2));
        if isfield(opt_set,'line_prob')
            if ischar(opt_set.line_prob)
                if strcmp(opt_set.line_prob,'all')
                    ndxL_prob = ones(length(mpc.branch(:,1)), size(cc,2));
                else
                    error('invalid entry for line_prob.');
                end
            else
                for rr = 1 : size(cc,2)
                    ndxL_prob(opt_set.line_prob(:,rr),rr) = 1;
                end
            end
        end
        if isfield(opt_set,'alpha')
            alpha = opt_set.alpha;
        else
            alpha = 0;
        end
        if isfield(opt_set,'correction')
            correction = opt_set.correction;
        else
            correction = 0;
        end
        if isfield(opt_set,'solver')
            solver = lower(opt_set.solver);
        else
            solver = 'null';
        end
        if isfield(opt_set,'tol_feas')
            tol_feas = opt_set.tol_feas;
        else
            tol_feas = 1e-6;
        end
        if isfield(opt_set,'tol_rank')
            tol_rank = opt_set.tol_rank;
        else
            tol_rank = 1e-6;
        end
    elseif nargin == 1
        ep = [0,0];
        cc{1} = double.empty(1,0);
        ndxL_prob = zeros(length(mpc.branch(:,1)), size(cc,2));
        alpha = 0;
        solver = 0;
        tol_feas = 1e-6;
        tol_rank = 1e-6;
        correction = 0;
    end
       
    nc = length(cc);

    [nb,busL,busN,busType,Pd,Qd,Vmax,Vmin,statusG,activeG,ng,genL,genN,incidentG,Qmax,Qmin,...
        Pmax,Pmin,activeL,statusL,nl,fbusL,tbusL,SlmMax,fbusN,tbusN,...
        incidentF,incidentT,Yf,Yt,YfP,YtP,Ybus,edges,c2,c1,c0] = Data_Reader(mpc,cc);

    [tw, perm, bags] = Permutation(busL,fbusL,tbusL,statusL,alpha);
    nBag = size(bags,1);
    
    if strcmp(solver,'null')
        if (tw <= 24) 
            solver = 'sdpt3';
        else
            solver = 'mosek';
        end
    end
    
    ndx = sparse(nb,nb);
    for kk = 1 : nBag
        ndx(busN(bags{kk}),busN(bags{kk})) = 1;
    end
    nlx = full(sum(sum(triu(ndx,1))));
    [fx,tx] = find(triu(ndx,1));


    disp('Initialization completed');

    if strcmp(solver,'mosek')
        MOSEK_Settings;
    end
    
    cvx_begin

    if strcmp(solver,'sdpt3')
        cvx_solver sdpt3
        cvx_solver_settings('maxit',120,'gaptol',1e-15,'inftol',1e-15,'steptol',1e-15);
    elseif strcmp(solver,'sedumi')
        cvx_solver sedumi
    elseif strcmp(solver,'mosek')
        cvx_solver mosek
    end

    cvx_precision best

    variable VV(nlx,nc) complex
    variable V2(nb,nc)
    variable Sg(ng,nc) complex

    expression sf(nl,nc)
    expression st(nl,nc)
    expression sfP(nl,nc)
    expression stP(nl,nc)
    expression Sb(nb,nc)
    expression W(nb*nc,nb*nc)

    for rr = 1 : nc
        in = (nb*(rr-1)+1) : (nb*rr);
        ex = setdiff(1 : nl, cc{rr});    

        W(in,in) = W(in,in) + sparse(fx, tx, VV(:,rr), nb, nb);
        W(in,in) = W(in,in) + sparse(tx, fx, conj(VV(:,rr)), nb, nb); 
        W(in,in) = W(in,in) + sparse(1:nb, 1:nb, V2(:,rr), nb, nb);

        sf(ex,rr) = conj(diag(Yf(ex, :) * W(in,in) * incidentF(ex, :).'));
        st(ex,rr) = conj(diag(Yt(ex, :) * W(in,in) * incidentT(ex, :).'));
        Sb(:,rr) = Pd + Qd*1i + conj(diag(Ybus{rr} * W(in,in)));

        sfP(ex,rr) = conj(diag(YfP(ex, :) * W(in,in) * incidentF(ex, :).'));
        stP(ex,rr) = conj(diag(YtP(ex, :) * W(in,in) * incidentT(ex, :).'));

        sf(cc{rr},rr) = zeros(length(cc{rr}),1);
        st(cc{rr},rr) = zeros(length(cc{rr}),1);

        sfP(cc{rr},rr) = zeros(length(cc{rr}),1);
        stP(cc{rr},rr) = zeros(length(cc{rr}),1);
    end

    obj = sum(c2.*real(Sg(:,1)).^2 + c1.*real(Sg(:,1)) + c0);
    panG = sum(sum(imag(Sg)));
    panQL = sum(sum(ndxL_prob .* abs(sfP + stP)));
    cost = obj + ep(1)*panG + ep(2)*panQL;
    
    minimize( cost );

    subject to

    for rr = 1 : nc
    for kk = 1 : nBag
        in = nb*(rr-1) + busN(bags{kk});
        W(in,in) == hermitian_semidefinite(size(bags{kk},2));
    end
    end

    Sb == incidentG.' * Sg;

    V2 <= (Vmax.^2) * ones(1,nc);
    V2 >= (Vmin.^2) * ones(1,nc);

    abs(sf) .* edges <= SlmMax * ones(1,nc);
    abs(st) .* edges <= SlmMax * ones(1,nc);

    real(Sg) <= Pmax * ones(1,nc);
    real(Sg) >= Pmin * ones(1,nc);
    imag(Sg) <= Qmax * ones(1,nc);
    imag(Sg) >= Qmin * ones(1,nc);

    for rr = 2 : nc
        abs(real(Sg(:,rr) - Sg(:,1))) <= correction;
    end

    cvx_end
    
    [updated_line_prob,bags_prob] = Rank_Check(W,bags,busN,fbusN,tbusN,nc,tol_rank);   
    
    [V_rec, Sg_rec, Sb_rec, sf_rec, st_rec, cost_rec] = Recovery(V2,W,Pd,Qd,ng, ...
        incidentG,Qmax,Qmin,Pmax,Pmin,correction,fbusN,tbusN,activeL,busType,Yf,Yt,Ybus,cc, ...
        solver,c2,c1,c0);
    
    disp('Cost value:');
    disp(cost_rec);
    disp('Number of high rank bags:');
    disp(length(bags_prob));
    disp('Number of problematic lines:');
    disp(length(updated_line_prob));
    
    [feas_flag, vio] = Feasibility_Check(V_rec,Sg_rec,Sb_rec,sf_rec,st_rec, ...
        incidentG,Vmax,Vmin,activeG,Qmax,Qmin,Pmax,Pmin,activeL,SlmMax,nc,tol_feas,correction);
    
    results.sdp.W = W;
    results.sdp.Wbag = cell(nc,1);
    for rr = 1 : nc
        in = (nb*(rr-1)+1) : (nb*rr);
        results.Wbag{rr} = W(in,in);
    end
    results.sdp.Sg = Sg;
    results.sdp.Sb = Sb;
    results.sdp.sf = sf;
    results.sdp.st = st;
    results.sdp.cost = obj;

    results.rec.V = V_rec;
    results.rec.Sg = Sg_rec;
    results.rec.Sb = Sb_rec;
    results.rec.sf = sf_rec;
    results.rec.st = st_rec;
    results.rec.cost = cost_rec;

    results.tw = tw;
    results.bags = bags;
    results.line_prob = updated_line_prob;
    results.bags_prob = bags_prob;
    results.violations = vio;
    results.feas_flag = feas_flag;

end
