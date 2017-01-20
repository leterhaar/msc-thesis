
clf
limits = [-10 10];
nul = [0 0];
plot3(limits, nul, nul, 'k');
hold on; grid on; box off;
plot3(nul, [0 20], nul, 'k');
plot3(nul, nul, limits, 'k');

size = [10 20 10];
origin = [0 0 0.02];
x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*size(1)+origin(1)+size(1)/2;
y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*size(2)+origin(2)+size(2)/2;
z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*size(3)+origin(3)+size(3)/2;
for i= [4 5]
    h=patch(x(:,i),y(:,i),z(:,i),'w');

    set(h,'LineStyle', 'none');
    set(h,'facecolor', 'g');
    set(h,'facealpha', 0.5);
end
% size = [2.5 2.5 2.5];
% origin = [0.0001 0 0];
% x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*size(1)+origin(1)+size(1)/2;
% y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*size(2)+origin(2)+size(2)/2;
% z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*size(3)+origin(3)+size(3)/2;
% for i=1:6
%     h=patch(x(:,i),y(:,i),z(:,i),'w');
%     hold on
%     set(h,'edgecolor','k')
%     set(h,'facecolor', 'w');
%     set(h,'facealpha', 0.8);
% end
axis('equal')
xlim(limits)
ylim([0 20])
zlim(limits)
a = 1.5;
b = 0.5;
n = 20;
s = linspace(0, 15, n);
% plot s
plot3(a+s*b, s, a+s*b, 'r')
plot3(a+s*b, s, zeros(n,1), 'g')
plot3(zeros(n,1), s, a+s*b, 'g')
% plot3(zeros(1,n), s, a+s*b)
set(gca, 'XTickLabels', [], 'YTickLabels', [], 'ZTickLabels', []);
text(12, 0, 0, 'x_1', 'fontsize', 22)
text(0, 20.5, 0.5, '\delta_i', 'fontsize', 22)
text(0, 0, 12, 'x_2', 'fontsize', 22)
text(a+s(end)*b, s(end), a+s(end)*b, 'a+\delta_ib', 'fontsize', 22, 'color', 'r');
view(42, 11)