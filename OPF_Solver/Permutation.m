% "Permutation.m" finds an elimination ordering of the network graph.

% SDP Solver of Optimal Power Flow: beta version
% Ramtin Madani(madani@ee.columbia.edu)
% Morteza Ashraphijuo(morteza.ashraphijuo@gmail.com)
% Javad Lavaei(lavaei@ee.columbia.edu)
% Columbia University
% Last Modified: October 07, 2014

function [tw, perm, bags] = Permutation(busL,fbusL,tbusL,statusL,alpha)

    nb = length(busL);
    nl = length(fbusL);
    
    busN = sparse(busL, ones(nb,1), (1:nb) );
    fbusN = busN(fbusL);
    tbusN = busN(tbusL);    

    perm = zeros(nb,1);
    neighbour = cell(nb,1);
    bags = cell(nb,1);

    for ii = 1 : nl
        if statusL(ii)
            neighbour{fbusN(ii)} = union( neighbour{fbusN(ii)}, tbusL(ii) );
            neighbour{tbusN(ii)} = union( neighbour{tbusN(ii)}, fbusL(ii) );
        end
    end

    fil = Fill_In(busN,neighbour,busL);
    deg = Deg(busN,neighbour,busL);

    tw = 1;
    for ii = 1 : nb
        [a,m] = min(fil);
        if a ~= 0
            [a,m] = min( alpha*deg + fil);
        end

        perm(ii) = busL(m);
        bb = union(perm(ii),neighbour{m});

        bags{ii} = bb;
        if  length(bb) > tw + 1
            tw = length(bb) - 1;
        end

        nn = busN(neighbour{m});
        nn2 = neighbour{m};
        for jj = 1 : length(nn)
            kk = nn(jj);
            nn2 = union(nn2,neighbour{kk});
            neighbour{kk} = setdiff(neighbour{kk},busL(m));
            neighbour{kk} = union(neighbour{kk},setdiff(neighbour{m},busL(kk)));
        end

        ff = Fill_In(busN,neighbour,nn2);
        fil(busN(nn2)) = ff;

        dd = Deg(busN,neighbour,neighbour{m});
        deg(nn) = dd;    

        fil(m) = Inf;
        deg(m) = Inf;

        neighbour{m} = 0;

    end
    disp('Upper bound on treewidth:');
    disp(tw);
end