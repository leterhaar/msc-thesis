function stretchY(ax, perc)
% stretches the y axis for more spacing
    ylims = ax.YLim;
    factor = perc/100;
    if ylims(1) > 0
        ax.YLim = [ylims(1)*(1-factor) ylims(2)*(1+factor)];
    else
        ax.YLim = [ylims(1)*(1+factor) ylims(2)*(1+factor)];
    end
end
