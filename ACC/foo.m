% make problem with P3 minimum set
t = 17;
Obj = AC_f(x, ac, wind, t);
C_det = AC_cons_det(x, ac, wind, t);

C_P3 = [];
for p = P3_act'
    i = p(1);
    j = p(2);
    
    C_P3 = [C_P3, AC_cons_scen(x, ac, wind.slice(i), t, j)];
end

opt = sdpsettings('solver', 'mosek', 'verbose', 0);
optimize([C_det, C_P3], Obj, opt);

AC_check(values_cell(x), ac, wind, t)

%%