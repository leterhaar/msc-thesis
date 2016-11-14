function h = initfig(name, id)
% for speed and brevity
% a function that clear a figure and sets it up for use
if nargin == 2
    h = figure(id);
else
    h = gcf;
end

if nargin >= 1
    h.Name = name;
end

clf
grid on
hold on
box on

h.Color = ones(3,1) * 250 / 255;
ax = gca;
ax.Color = ones(3,1) * 245 / 255;
% ax.FontName = 'Fira Sans';
% ax.XLabel.FontSize = 18;
% ax.YLabel.FontSize = 18;


assignin('base', 'blue', [36 55 58]/255);
assignin('base', 'green', [128 205 193]/255); 
assignin('base', 'orange', [235 128 27]/255);
assignin('base', 'brown', [191 128 64]/255);	
end