function psd_cone
addpath('../misc');
initfig('PSD Cone', 1);
hold off
handle = plotCone([0 0 0], ones(3,1)*0.8);
A = [0.5 0.4 0.5];
verify(is_psd([A(1) A(2); A(2) A(3)]), 'A not PSD');
pointA = plot3(A(1), A(2), A(3), ...
               'ro', 'markersize', 5,'markerfacecolor', 'r');
handle2 = plotCone(A, 'r');
legend([handle, handle2], 'PSD cone', 'PSD cone for A+\delta B');
xlim([0 1]);
ylim([-1 1]);
zlim([0 1]);
xlabel('x');
ylabel('y');
zlabel('z');
view(-45, 20);
end

function handle = plotCone(orig, color)
step = 0.1;
eps = 0.01;
for y = step:step:2
    for x = step:step:2
        z = y^2/x;
        nextz = y^2 / (x+step);
        prevz = y^2 / (x-step);
        handle = patch('xdata', [0 x x+step 0]+orig(1), ...
      'ydata', [0 y y 0]+orig(2), ...
      'zdata', [0 z nextz 0]+orig(3), 'edgealpha', 0, ...
      'facecolor', color, 'facealpha', 0.8);
        patch('xdata', [0 x x+step 0]+orig(1), ...
      'ydata', [0 -y -y 0]+orig(2), ...
      'zdata', [0 z nextz 0]+orig(3), 'edgealpha', 0, ...
      'facecolor', color, 'facealpha', 0.8);
        if rem(x,.3) == 0 && rem(y, .3) == 0
            plot3(orig(1)+[0 x], orig(2)+[0 y], orig(3)+[0 z], 'k', 'linewidth', .1);
            plot3(orig(1)+[0 x], orig(2)+[0 -y], orig(3)+[0 z], 'k', 'linewidth', .1);
        end
        hold on
    end
end
x = -1:0.1:1;
plot3(-orig(1) + x.^2, -orig(2)+x, orig(3)+ones(length(x),1), 'k', 'linewidth', .8);
plot3(-orig(1) + -ones(length(x),1), orig(2)+ x, orig(3) + x.^2, 'k', 'linewidth', .8);
end
