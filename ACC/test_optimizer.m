% MWE for weird optimizer behaviour
N = 10;
d = 5;
deltas = rand(N,d);
x = sdpvar(d,1);
delta = sdpvar(1,1);

