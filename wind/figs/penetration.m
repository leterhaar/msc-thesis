% read in data
data = csvread('gwec.csv')';
% convert to gigawatt
data = data ./ 1000;
annual = data(:,1);
cumulative = data(:,2);

% create tick labels
years = 0:15;
labels = cell(16,1);
for year = years
    if rem(year, 2) == 0

        if year < 10
            labels{year+1} = sprintf('`0%i', year);
        else
            labels{year+1} = sprintf('`%i', year);
        end
    else
        labels{year+1} = '';
    end
    
end


% check if cumulative makes sense
cum_check = zeros(length(annual),1);
for i = 1:16
    cum_check(i) = sum(annual(1:i));
end
cum_check = cum_check + cumulative(1);


addpath('../../misc');
initfig('wind-penetration',1);
yyaxis left
b = bar(years, cumulative);
b.FaceColor = blue;
ylabel('Cumulative installed wind capacity [GW]');
xlim([-0.75 15.75]);
set(gca, 'XTick', years, 'XTickLabels', labels, 'YColor', blue);

yyaxis right
l = plot(years, annual, 'o-', 'color', green);
l.MarkerFaceColor = l.Color;
set(gca, 'YColor', green);
ylabel('Annual installed wind capacity [GW]'); 

title('Global installed wind capacity 2000-2015')