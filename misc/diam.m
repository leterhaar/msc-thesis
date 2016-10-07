function [dm, d] = diam(A)
% find the diameter of a graph corresponding to adjacancy matrix A
G = digraph(A);
d = distances(G);
assert(all(all(isfinite(d))), 'Graph not fully connected!');
dm = max(d(:));