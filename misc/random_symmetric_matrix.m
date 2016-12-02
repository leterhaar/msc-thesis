function M = random_symmetric_matrix(N)
% returns a random symmetric matrix of dimensions NxN
% from https://nl.mathworks.com/matlabcentral/answers/123643-how-to-create-a-symmetric-random-matrix
    d = 1000000*rand(N,1); % The diagonal values
    t = triu(bsxfun(@min,d,d.').*rand(N),1); % The upper trianglar random values
    M = diag(d)+t+t.'; % Put them together in a symmetric matrix
end
