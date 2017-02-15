function Xres = project_PSD(X, sparsity)
% project a matrix in the PSD cone through eigenvector decompositiobn

    if nargin < 2
        sparsity = X ~= 0;
    end
    d = size(X,1);
    if d == 1
        Xres = max(real(X), 0);
        return
    end
    
    % shrink it down to a small matrix
    nzs = find(sum(sparsity)); % find the non-zero columns / rows
    nnzs = length(nzs);
    selector = zeros(d, nnzs);
    for i = 1:nnzs
        selector(nzs(i), i) = 1;
    end

    X_small = selector' * X * selector;
    try
        % do an eigenvalue decomposition
        [Q, D] = eig(X_small);
        % remove negative eigenvalues
        D = max(real(D), 0);          
        % reconstruct X_small
        X_small = Q * D / Q;
    catch e
        e.getReport()
        keyboard
    end
    Xres = selector * X_small * selector';
end