function divided = sparse_division(D,C,E)
% see "ADMM for Sparse Semidefinite Programming with Applications to
% Optimal Power Flow Problem" by Madani et al. - Notation 1

    divided = D ./ E; % entrywise division
    C(C == 0) = nan;
    divided = divided ./ C;
    divided(isnan(divided)) = 0;
end