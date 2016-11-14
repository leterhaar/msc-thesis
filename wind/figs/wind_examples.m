% cd to one dir up
% cd ../
wind = wind_model(dc, 24, 0.2);
wind.generate(150);

initfig('Wind', 2);

h = plot(wind.P_w, 'linewidth', 2, 'color', blue);
xlim([1 24]);
for line = h'
    line.Color = line.Color + ones(1,3)*rand()/2;
end

xlabel('hours');
ylabel('P_w');
title('Historic wind realizations');