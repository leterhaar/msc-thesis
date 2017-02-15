function C = PSD_bags(W, bags)
% defines PSD constraints on the entries of W defined by bags
    C = [];
    for i = 1:length(bags)
        bag = bags{i};
        C = [C, W(bag, bag) >= 0];
    end
end