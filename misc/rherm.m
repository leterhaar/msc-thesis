function [ A ] = rherm(n)
% Generate a hermitian matrix with random values.
% rherm(n) returns a random hermitian matrix of size n*n
if exist('n')~=1, n=2; end
A = [zeros(1,0),(2*rand(1,1)-1),(2*rand(1,n-1)-1)+(2*rand(1,n-1)-1)*1i;zeros(n-1,n)];
for i=1:n-1, A = A + [zeros(n-i,n);zeros(1,n-i),(2*rand(1,1)-1),(2*rand(1,i-1)-1)+(2*rand(1,i-1)-1)*1i;zeros(i-1,n)];end
for i=1:n-1, for j=1:n-i, A = A + [zeros(n-i,n);zeros(1,j-1),conj(A(j,n-i+1)),zeros(1,n-j);zeros(i-1,n)];end;end
end