clf
N = 24;
t = linspace(0, 24, N);
y = 0.5 + 0.5 * sin(t*pi/24) + 0.05*randn(1,N);

bar(t,y);
xlim([0.5,24.5])
