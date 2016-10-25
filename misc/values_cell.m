function xstar = values_cell(x)
% maps the value function to the elements of x:
%
% values_cell(x) = {value(x{1}), value(x{2}), ... value(x{n})}
%
% useful for cells of sdpvars

n = length(x);
xstar = cell(n,1);
for i = 1:n
    xstar{i} = value(x{i});
end
