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
% ax.Color = ones(3,1) * 245 / 255;
% ax.FontName = 'Fira Sans';
% ax.XLabel.FontSize = 18;
% ax.YLabel.FontSize = 18;
set(h, 'defaultLineLineWidth', 2)
title(name);

% Presentation colors
% assignin('caller', 'blue', [36 55 58]/255);
% assignin('caller', 'green', [128 205 193]/255); 
% assignin('caller', 'orange', [235 128 27]/255);
% assignin('caller', 'brown', [191 128 64]/255);	
assignin('caller', 'blue', [0    0.4470    0.7410] );
assignin('caller', 'green', [ 0.8500    0.3250    0.0980]); 
assignin('caller', 'orange', [0.9290    0.6940    0.1250] );
assignin('caller', 'brown', [ 0.4940    0.1840    0.5560]);	

end