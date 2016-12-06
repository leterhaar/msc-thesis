% a fake process that takes some time
p = progress('Testing', 7200);
for i = 1:2000
    pause(1);
    p.ping();
end