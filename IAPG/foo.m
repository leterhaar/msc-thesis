clc
N = 4;
all_is = 1:N;

for k = N+1:20
    i = rem(k-1, N) + 1;
    for j = all_is(all_is ~= i)
        l = floor(k/N)*N - (N-j);
        fprintf('k%i i%i j%i l%i\n',k, i, j, l);
    end
end

