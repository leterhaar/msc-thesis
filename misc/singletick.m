ax = gca;
xlims = xlim;

% get separation between ticks
sep = ax.XTick(2) - ax.XTick(1);
if sep < 1
    ax.XTick = xlims(1):xlims(2);
end