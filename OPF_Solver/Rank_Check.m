% "Rank_Check.m" checks the rank of each bag and outputs problematic bags/lines.

% SDP Solver of Optimal Power Flow: beta version
% Ramtin Madani(madani@ee.columbia.edu)
% Morteza Ashraphijuo(morteza.ashraphijuo@gmail.com)
% Javad Lavaei(lavaei@ee.columbia.edu)
% Columbia University
% Last Modified: October 07, 2014

function [line_prob,bags_prob] = Rank_Check(W,bags,busN,fbusN,tbusN,nc,tol_rank)

    nb = size(W,1)/nc;
    nl = length(fbusN);
    nBag = size(bags,1);
    
    ndxBags_prob = zeros(nBag,nc);
    ndx_prob = sparse(nc*nb, nc*nb);
    for rr = 1 : nc
    for kk = 1 : nBag
        if length(bags{kk})>1
            in = nb*(rr-1);
            nodes = busN(bags{kk});
            ee = eigs(W(in + nodes, in + nodes));
            if ee(2) > tol_rank
                ndxBags_prob(kk,rr) = 1;
                ndx_prob(in + nodes, in + nodes) = 1;
            end
        end
    end
    end

    [fx_prob,tx_prob] = find(triu(ndx_prob,1));
    ndxL_prob = zeros(nl,nc);
    for rr = 1 : nc
    for ii = 1 : nl
        in = nb*(rr-1);
        uu = intersect(find(fx_prob == (fbusN(ii) + in)), find(tx_prob == (tbusN(ii) + in)));
        dd = intersect(find(fx_prob == (tbusN(ii) + in)), find(tx_prob == (fbusN(ii) + in)));
        if (~isempty(uu) || ~isempty(dd))
            ndxL_prob(ii, rr) = 1;
        end
    end
    end
    line_prob = find(ndxL_prob);
    bags_prob = find(ndxBags_prob);
    
end