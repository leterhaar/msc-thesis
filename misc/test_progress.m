% a fake process that takes some time
p = progress('Testing', 100);
for i = 1:100
    pause(0.01);
    p.ping();
end