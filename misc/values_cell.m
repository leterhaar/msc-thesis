function xstar = values_cell(x)
% maps the value function to the elements of x:
%
% values_cell(x) = {value(x{1}), value(x{2}), ... value(x{n})}
%
% useful for cells of sdpvars

[n,m] = size(x);
xstar = cell(n,m);
for i = 1:n
    for j = 1:m
        xstar{i,j} = value(x{i,j});
    end
end
