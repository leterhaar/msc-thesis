function assign_cell(x, values)
% maps the assign function to two cell arrays
n = length(x);
for i = 1:n
    assign(x{i}, values{i});
end