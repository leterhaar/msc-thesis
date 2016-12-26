function res = cell_divide(X, a)
% returns every element of the cell elemenet wise divided by a
    res = cell(length(X), 1);
    for i = 1:length(X)
        res{i} = X{i} ./ a;
    end
end