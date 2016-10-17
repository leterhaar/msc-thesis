function adj = random_graph(N, dm, method)
% heuristic script that creates a random graph with N vertices 
% with diameter diam

% debug
if nargin == 0
    N = 30;
    dm = 4;
end

if nargin < 3
    method = 'rgg';
end

assert(N > 2, 'Too small N');

adj = zeros(N);

if strcmpi(method, 'rand')
    % to create a fully connected graph, every node should at least connect to
    % one other connected node
    % connect the first node to some other node
    to = randi(N-1)+1;
    adj(1, to) = 1;
    adj(to, 1) = 1;
    connected_nodes = [1, to, NaN(1, N-1)];

    for i = 2:N

        % pick randomly from the connected nodes
        while 1
            to = connected_nodes(randi(i));
            if to ~= i % cannot connect to self
                break
            end
        end

        connected_nodes(i+1) = i;
        adj(i, to) = 1;
        adj(to, i) = 1;
    end

    % now a fully connected graph is created
    % keep adding edges until diameter is reached
    while diam(adj) > dm
        adj(randi(N), randi(N)) = 1;
    end
    
    % remove self loops
    adj = adj - diag(diag(adj));

elseif strcmpi(method, 'rgg')
    
    solution_found = 0;
    
    while not(solution_found)

        % create N points on a 2D plane
        points = rand(N,2);
        dist = zeros(N);

        % calculate distance between points
        for i = 1:N
            diff = points-(ones(N,1)*points(i, :));
            dist(:, i) = sqrt(diff(:,1).^2 + diff(:,2).^2);
        end

        dist = dist + diag(inf(N,1));
        thresholds = 1./(1:0.1:100);

        try
            for threshold = thresholds

                adj = dist <= threshold;
                adj = adj - diag(diag(adj));
                if diam(adj) == dm
                    solution_found = 1;
                    break
                end
            end
        catch
            solution_found = 0;
        end
    end
end