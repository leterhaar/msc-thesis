function res = all_close(A, B, tol)
% compares A and B with tolerance

    if nargin < 3
        tol = 1e-6;
    end

    assert(all(size(A) == size(B)), 'Should be of same dimensions');

    % loop over cell with doubles
    if isa(A, 'cell') && isa(B, 'cell')

        res = 1;
        for i = 1:length(A)
            A_tmp = A{i};
            B_tmp = B{i};
            
            assert(all(size(A_tmp) == size(B_tmp)), 'Should be of same dimensions');
            
            % vectorize and check difference
            if any(abs(A_tmp(:) - B_tmp(:)) >= tol)
                res = 0;
            end
        end
    
    else
        % vectorize and check difference
        res = all(abs(A(:) - B(:)) < tol);
    end
    
    res = logical(res);
end