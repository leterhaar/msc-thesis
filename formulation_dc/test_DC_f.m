% create objective function with sdpvar

Obj_test = DC_f(x, dc, wind);

assert(is_scalar(Obj_test), 'wrong dims');
assert(isa(Obj_test, 'sdpvar'), 'wrong type');

% create objective function with doubles

Obj_test = DC_f(random_x, dc, wind);

assert(is_scalar(Obj_test), 'wrong dims');
assert(isa(Obj_test, 'double'), 'wrong type');
assert(not(isnan(Obj_test)), 'has nans');

% check if this is the same as summing all individual objective functions
% as found in previous notation
Obj_old = 0;


for t = 1:N_t
    P_G = random_x(1:dc.N_G, t);
    R_us = random_x(dc.N_G+1:2*dc.N_G, t);
    R_ds = random_x(2*dc.N_G+1:3*dc.N_G, t);
    % loop over generators and add generator cost at time t
    for k = 1:dc.N_G
        Obj_old = Obj_old + dc.c_qu(k) * (P_G(k))^2 + ...
                                dc.c_li(k) * P_G(k);
    %       Obj = Obj + dc.c_us(k) * P_G(k);
    end

    % add reserve requirements costs
    Obj_old = Obj_old + (dc.c_us' * R_us + dc.c_ds' * R_ds);
end

assert(all_close(Obj_old, Obj_test, 1e-5), 'Not the same');