function res = cell_add(A, B)
% adds the elements of two cells
    assert(length(A) == length(B), 'Should be same length');
    res = cell(length(A), 1);
    for i = 1:length(A)
        assert(all(size(A) == size(B)), 'Should have same elements');
        res{i} = A{i} + B{i};
    end
end