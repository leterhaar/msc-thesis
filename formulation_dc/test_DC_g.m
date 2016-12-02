dim_g = (4*dc.N_G + 2*dc.N_l)*N_t;

%% test with sdp vars and no j_des
g = DC_g(x, dc, wind.slice(1));

% should return an sdpvar of the right dimensions
assert(isa(g, 'sdpvar'), 'g is not an sdp var');
assert(all(size(g) == [dim_g, 1]), 'g is not the right dimension');

%% test with sdpvars and a j_des
random_j = nan(1,10);
for i = 1:10
    random_j(i) = randi(dim_g);
end

for j = random_j
    g = DC_g(x, dc, wind.slice(1), j); 
    
    % should return an sdpvar of the right dimensions
    assert(isa(g, 'sdpvar'), 'g is not an sdp var');
    assert(is_scalar(g), 'g is not the right dimension');
end



%% test with actual values and no j_des
g = DC_g(random_x, dc, wind.slice(1));

% should return a double 
assert(isa(g, 'double'), 'g is not a double');
% without nans
assert(not(any(isnan(g))), 'g has nans');
% of the right size
assert(all(size(g) == [dim_g, 1]), 'g is not the right dimension');

%% test with random values and j_des
for j = random_j
    % retrieve the residual
    g = DC_g(random_x, dc, wind.slice(1), j);

    % check dimensions and types
    assert(isa(g, 'double'), 'g is not a double');
    assert(all(size(g) == [1, 1]), 'g is not the right dimension');
    assert(not(isnan(g)), 'Should not be nan');
end