Ztilde = rherm(10);
while is_psd(Ztilde)
    Ztilde = rherm(10);
end

X = sdpvar(10, 10, 'hermitian', 'complex');
optimize([X >= 0], norm(X-Ztilde, 'fro')^2, sdpsettings('verbose', 0));

Zopt = value(X);
Zopt2 = project_PSD(Ztilde);
assert(all_close(Zopt, Zopt2, 1e-4))


