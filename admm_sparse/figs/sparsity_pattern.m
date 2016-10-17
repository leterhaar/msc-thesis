clf
addpath('../../misc');

assert(exist('sparsity', 'var') == 1, 'Make sure sparsity is in workspace first!');

subplot(121)
spy(sparsity)
title('Sparsity pattern');

ax = subplot(122);
plot(graph(sparsity))
title('Graph of sparsity pattern');
ax.XTick = [];
ax.YTick = [];
axis('square');