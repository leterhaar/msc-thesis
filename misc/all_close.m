function [close, infnorm] = all_close(A, B, tol)
% [res, infnorm] = all_close(A, B, tol)
%
% compares A and B with absolute tolerance tol
%
% RETURNS
% =======
% close: logical to indicate whether A and B are close
% infnorm: infity norm (max abs) of difference between A and B 

    if nargin < 3
        tol = 1e-6;
    end

    assert(all(size(A) == size(B)), 'Should be of same dimensions');
    close = 1;
    
    % loop over cell with doubles
    if isa(A, 'cell') && isa(B, 'cell')

        infnorm = -inf;
        for i = 1:length(A)
            A_tmp = A{i};
            B_tmp = B{i};
            
            assert(all(size(A_tmp) == size(B_tmp)), 'Should be of same dimensions');
            
            % vectorize and check difference
            if any(abs(A_tmp(:) - B_tmp(:)) >= tol)
                close = 0;
            end
            
            infnorm_tmp = norm(A_tmp(:) - B_tmp(:), 'inf');
        end
    
        % take the maximum
        infnorm = max(infnorm_tmp, infnorm);
    else
        % vectorize and check difference
        if any(abs(A(:) - B(:)) >= tol)
            close = 0;
        end
        infnorm = norm(A(:) - B(:), 'inf');
    end
    
    close = logical(close);
    
    
end