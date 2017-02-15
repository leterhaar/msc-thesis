% "Deg.m" computes degree of the input nodes in the network graph.

% SDP Solver of Optimal Power Flow: beta version
% Ramtin Madani(madani@ee.columbia.edu)
% Morteza Ashraphijuo(morteza.ashraphijuo@gmail.com)
% Javad Lavaei(lavaei@ee.columbia.edu)
% Columbia University
% Last Modified: October 07, 2014

function [ff] = Deg(busN,neighbour,bb)
    ff = zeros(length(bb),1);
    for ii = 1 : length(bb)
        qq = busN(bb(ii));
        ff(ii) = length(neighbour{qq});
    end
end