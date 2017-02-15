function [Wopt, Vopt] = completeW(W)
% extracts the voltage profile from an incomplete W and reconstructs W
    
    % extract network information from caller
    ac = evalin('caller', 'ac');
    for bag = ac.bags'
        if svd_rank(W(bag{:},bag{:}), 1e-1) > 1
            fprintf('Submatrix with buses %s not rank 1', num2str(bag{:}));
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