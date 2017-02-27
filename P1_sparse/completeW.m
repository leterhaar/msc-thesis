function [Wopt, Vopt] = completeW(W, ac)
% extracts the voltage profile from an incomplete W and reconstructs W
% 
% ARGS
% =====
% W : (partially filled) state matrix
% ac : instance of network model (optional) 
%      if left empty, 'ac' will be evaluated in caller
    
    % extract network information from caller
    if nargin < 2
        ac = evalin('caller', 'ac');
    end
    for bag = ac.bags'
        the_rank = svd_rank(W(bag{:},bag{:}));
%         svds(W(bag{:},bag{:}))
        if the_rank > 1
            fprintf('Submatrix with buses %s is rank %i!\n', num2str(bag{:}), the_rank);
        end
    end
    
    % get diagonal
    X = sqrt(diag(W));
    Vmag = X(1:ac.N_b);
    Vang = X(1+ac.N_b:end);
    
    
    % find sign for every bus
    signs = nan(ac.N_b, 1);
    
    % make connection matrix
    lines = zeros(ac.N_b);
    for l = 1:ac.N_l
        i = ac.from_to(l,1);
        j = ac.from_to(l,2);
        lines(i,j) = 1;
        lines(j,i) = 1;
    end
    
    
    for to_bus = 1:ac.N_b
        from_bus = find(lines(to_bus, :),1);
        signs(to_bus) = sign(W(from_bus, ac.N_b+to_bus));
    end
    
    X = [Vmag; signs.*Vang];
    Wopt = X * X';
    Vopt = Vmag + sqrt(-1)*(signs.*Vang);
end