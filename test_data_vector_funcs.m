% Test make_data_vectors.m and unpack_data_vectors.m

%% Test make_data_vectors

vdata = make_data_vectors({0});
assert(isequal(vdata.vecs, {0}))
assert(isequal(vdata.types, {'double'}))
assert(isequal(vdata.dims, {[1 1]}))
assert(isequal(vdata.n_els, {[1]}))

vdata = make_data_vectors({1, 2, 3});
assert(isequal(vdata.vecs, {1, 2, 3}))
assert(isequal(vdata.types, {'double', 'double', 'double'}))
assert(isequal(vdata.dims, {[1 1], [1 1], [1 1]}))
assert(isequal(vdata.n_els, {1, 1, 1}))

vdata = make_data_vectors({1, [2; 3], [4 6 8; 5 7 9]});
assert(isequal(vdata.vecs, {1, [2 3], [4 5 6 7 8 9]}))
assert(isequal(vdata.types, {'double', 'double', 'double'}))
assert(isequal(vdata.dims, {[1 1], [2 1], [2 3]}))
assert(isequal(vdata.n_els, {1, 2, 6}))

x1 = {1 3 5; 2 4 6};
vdata_x1 = make_data_vectors({x1});
assert(isequal(vdata_x1.vecs, {1:6}))
assert(isequal(vdata_x1.types, {{'double', 'double', 'double'; 'double', 'double', 'double'}}))
assert(isequal(vdata_x1.dims, {{[1 1], [1 1], [1 1]; [1 1], [1 1], [1 1]}}))
assert(isequal(vdata_x1.n_els, {6}))

x2 = {[7 8 9], 10, [11; 12], [13 15 17; 14 16 18]};
vdata_x2 = make_data_vectors({x2});
assert(isequal(vdata_x2.vecs, {7:18}))
assert(isequal(vdata_x2.types, {{'double', 'double', 'double', 'double'}}))
assert(isequal(vdata_x2.dims, {{[1 3], [1 1], [2 1], [2 3]}}))
assert(isequal(vdata_x2.n_els, {12}))

vdata_y = make_data_vectors({x1, x2, [19; 20; 21]});
assert(isequal(cell2mat(vdata_y.vecs), 1:21))
assert(isequal(vdata_y.types, {vdata_x1.types{1}, vdata_x2.types{1}, 'double'}))
assert(isequal(vdata_y.dims, {vdata_x1.dims{1}, vdata_x2.dims{1}, [3 1]}))
assert(isequal(vdata_y.n_els, {6, 12, 3}))

A = [1 3 5; 2 4 6];
vdata = make_data_vectors({A, x2, 19, [20; 21]});
assert(isequal(cell2mat(vdata.vecs), 1:21))
assert(isequal(vdata.types, {'double', vdata_x2.types{1}, 'double', 'double'}))
assert(isequal(vdata.dims, {size(A), vdata_x2.dims{1}, [1 1], [2 1]}))
assert(isequal(vdata.n_els, {6, 12, 1, 2}))


%% Test unpack_data_vectors

vdata = make_data_vectors({1});
vars = unpack_data_vectors(vdata);
assert(isequal(vars, {1}))

vdata = make_data_vectors({1, 2});
vars = unpack_data_vectors(vdata);
assert(isequal(vars, {1, 2}))

vdata = make_data_vectors({1, [2; 3], [4 5]});
vars = unpack_data_vectors(vdata);
assert(isequal(vars, {1, [2; 3], [4 5]}))

% Cell array
a = 1;
b = [2 3; 4 5];
c = {6, [7; 8], 9};
vdata = make_data_vectors({a, b, c});
vars = unpack_data_vectors(vdata);
assert(isequal(vars, {a, b, c}))

% Cell array
a = 1;
b = [2 3; 4 5];
c = {6, [7; 8], 9};
d = {10; 11};
vdata = make_data_vectors({a b; c d});
vars = unpack_data_vectors(vdata);
assert(isequal(vars, {a b; c d}))

% Nested cell arrays
a = 1;
b = [2 3; 4 5];
c = {6, {[7 8 9], [10; 11]}, 12};
vdata = make_data_vectors({a, b, c});
vars = unpack_data_vectors(vdata);
assert(isequal(vars, {a, b, c}))


%% Test with different numeric data types

i = int16(10);
vdata = make_data_vectors({i}, 'int16');
assert(isequal(vdata.vecs, {int16(10)}))
assert(isequal(vdata.types, {'int16'}))
assert(isequal(vdata.dims, {[1 1]}))

vdata = make_data_vectors({uint8(1), uint8([2; 3]), ...
    uint8([4 6 8; 5 7 9])}, 'uint8');
assert(isequal(vdata.vecs, {uint8(1), uint8([2 3]), uint8([4 5 6 7 8 9])}))
assert(isequal(vdata.types, {'uint8', 'uint8', 'uint8'}))
assert(isequal(vdata.dims, {[1 1], [2 1], [2 3]}))


%% Test with Kalman filter object

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Sequenece length
nT = 20;

% Simulate system
X0 = zeros(n,1);
t = Ts*(0:nT)';
Wp = sample_random_shocks(nT+1, epsilon, sigma_wp(2), sigma_wp(1));
U = zeros(nT+1,1);
U(t>=5) = 1;
Y = lsim(Gpss,[U Wp],t,X0);
V = sigma_M*randn(nT+1, 1);
Ym = Y + V;  % measurement

% Make copies of KF1
obs1 = KF1.copy();
obs2 = KF1.copy();

% Convert dynamic variables to vdata struct
vars = get_obs_vars(obs2);
vdata = make_data_vectors(struct2cell(vars)');

% Simulate observer
Xk_est = {nan(nT+1,n), nan(nT+1,n)};
Yk_est = {nan(nT+1,ny), nan(nT+1,ny)};
for i = 1:nT

    uk = U(i,:)';
    yk = Ym(i,:)';
    obs1.update(yk, uk);

    % Save observer estimates
    Xk_est{1}(i+1,:) = obs1.xk_est';
    Yk_est{1}(i+1,:) = obs1.yk_est';

    % Unpack vdata struct and reconstruct observer from KF1
    vars = unpack_data_vectors(vdata);
    assert(isequal(vars{1}, obs2.xkp1_est))
    assert(isequal(vars{2}, obs2.Pkp1))
    obs2 = KF1;  % makes a new copy
    obs2.xkp1_est = vars{1};
    obs2.Pkp1 = vars{2};

    % Observer updates
    obs2.update(yk, uk);

    % Convert dynamic variables to vdata struct
    vars = get_obs_vars(obs2);
    vdata = make_data_vectors(struct2cell(vars)');

    % Save observer estimates
    Xk_est{2}(i+1,:) = obs2.xk_est';
    Yk_est{2}(i+1,:) = obs2.yk_est';

end

% Check observer estimates are identical
assert(isequaln(Xk_est{1}, Xk_est{2}))
assert(isequaln(Yk_est{1}, Yk_est{2}))


%% Test with MKF observer object

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Sequenece length
nT = 20;

% Simulate system
X0 = zeros(n,1);
t = Ts*(0:nT)';
Wp = sample_random_shocks(nT+1, epsilon, sigma_wp(2), sigma_wp(1));
U = zeros(nT+1,1);
U(t>=5) = 1;
Y = lsim(Gpss,[U Wp],t,X0);
V = sigma_M*randn(nT+1, 1);
Ym = Y + V;  % measurement

% Make copies of MKF_SP1
obs1 = MKF_SP1.copy();
obs2 = MKF_SP1.copy();

% Convert dynamic variables to vdata struct
vars = get_obs_vars(obs2);
vars_double = {vars.xkp1_est, vars.Pkp1, vars.p_seq_g_Yk, ...
    vars.xkp1_est_f, vars.Pkp1_f};
vdata = make_data_vectors(vars_double);
vdata_int16 = make_data_vectors(struct2cell(vars.int16)', 'int16');

% Simulate observer
Xk_est = {nan(nT+1,n), nan(nT+1,n)};
Yk_est = {nan(nT+1,ny), nan(nT+1,ny)};
for i = 1:nT

    uk = U(i,:)';
    yk = Ym(i,:)';
    obs1.update(yk, uk);

    % Save observer estimates
    Xk_est{1}(i,:) = obs1.xk_est';
    Yk_est{1}(i,:) = obs1.yk_est';

    % Unpack vdata struct and copy values back to observer
    vars_double = unpack_data_vectors(vdata);
    vars_int16 = unpack_data_vectors(vdata_int16);
    vars = struct();
    vars.xkp1_est = vars_double{1};
    vars.Pkp1 = vars_double{2};
    vars.p_seq_g_Yk = vars_double{3};
    vars.xkp1_est_f = vars_double{4};
    vars.Pkp1_f = vars_double{5};
    vars.int16.rk = vars_int16{1};
    vars.int16.f_main = vars_int16{2};
    vars.int16.f_hold = vars_int16{3};

    % Build a copy from the variable data
    obs2_restored = set_obs_vars(MKF_SP1.copy(), vars);
    assert(isequal(obs2_restored.xkp1_est, obs2.xkp1_est))
    assert(isequal(obs2_restored.Pkp1, obs2.Pkp1))
    assert(isequal(obs2_restored.p_seq_g_Yk, obs2.p_seq_g_Yk))
    assert(isequal(obs2_restored.filters.Xkp1_est, obs2.filters.Xkp1_est))
    assert(isequal(obs2_restored.filters.Pkp1, obs2.filters.Pkp1))
    assert(isequal(obs2_restored.rk, obs2.rk))
    assert(isequal(obs2_restored.f_main, obs2.f_main))
    assert(isequal(obs2_restored.f_hold, obs2.f_hold))

    % Observer updates
    obs2.update(yk, uk);

    % Re-convert dynamic variables to vdata struct for next
    % time step
    vars = get_obs_vars(obs2);
    vars_double = {vars.xkp1_est, vars.Pkp1, vars.p_seq_g_Yk, ...
        vars.xkp1_est_f, vars.Pkp1_f};
    vdata = make_data_vectors(vars_double);
    vdata_int16 = make_data_vectors(struct2cell(vars.int16)', 'int16');

    % Save observer estimates
    Xk_est{2}(i,:) = obs2.xk_est';
    Yk_est{2}(i,:) = obs2.yk_est';

end

% Check observer estimates are identical
assert(isequaln(Xk_est{1}, Xk_est{2}))
assert(isequaln(Yk_est{1}, Yk_est{2}))


%% Test with MKF_SF observer object

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Sequenece length
nT = 20;

% Simulate system
X0 = zeros(n,1);
t = Ts*(0:nT)';
Wp = sample_random_shocks(nT+1, epsilon, sigma_wp(2), sigma_wp(1));
U = zeros(nT+1,1);
U(t>=5) = 1;
[Y,T,X] = lsim(Gpss,[U Wp],t,X0);
V = sigma_M*randn(nT+1, 1);
Ym = Y + V;  % measurement

% Make copies of MKF_SP1
obs1 = MKF_SF1.copy();
obs2 = MKF_SF1.copy();

% Convert dynamic variables to vdata struct
vars = get_obs_vars(obs2);
vars_double = {vars.xkp1_est, vars.Pkp1, vars.p_seq_g_Yk, ...
    vars.xkp1_est_f, vars.Pkp1_f};
vdata = make_data_vectors(vars_double);
vdata_int16 = make_data_vectors(struct2cell(vars.int16), 'int16');

% Simulate observer
Xk_est = {nan(nT+1,n), nan(nT+1,n)};
Yk_est = {nan(nT+1,ny), nan(nT+1,ny)};
for i = 1:nT

    uk = U(i,:)';
    yk = Ym(i,:)';
    obs1.update(yk, uk);

    % Save observer estimates
    Xk_est{1}(i+1,:) = obs1.xk_est';
    Yk_est{1}(i+1,:) = obs1.yk_est';

    % Unpack vdata struct and copy values back to observer
    vars_double = unpack_data_vectors(vdata);
    vars_int16 = unpack_data_vectors(vdata_int16);
    vars = struct();
    vars.xkp1_est = vars_double{1};
    vars.Pkp1 = vars_double{2};
    vars.p_seq_g_Yk = vars_double{3};
    vars.xkp1_est_f = vars_double{4};
    vars.Pkp1_f = vars_double{5};
    vars.int16.rk = vars_int16{1};
    vars.int16.i = vars_int16{2};
    vars.int16.i_next = vars_int16{3};
    vars.int16.i2 = vars_int16{4};
    vars.int16.i2_next = vars_int16{5};
    vars.int16.seq = vars_int16{6};

    % Build a copy from the variable data

    obs2_restored = set_obs_vars(MKF_SF1.copy(), vars);  % makes a new copy
    assert(isequal(obs2_restored.xkp1_est, obs2.xkp1_est))
    assert(isequal(obs2_restored.Pkp1, obs2.Pkp1))
    assert(isequal(obs2_restored.p_seq_g_Yk, obs2.p_seq_g_Yk))
    assert(isequal(obs2_restored.filters.Xkp1_est, obs2.filters.Xkp1_est))
    assert(isequal(obs2_restored.filters.Pkp1, obs2.filters.Pkp1))
    assert(isequal(obs2_restored.rk, obs2.rk))
    assert(isequal(obs2_restored.i, obs2.i))
    assert(isequal(obs2_restored.i_next, obs2.i_next))

    % Observer updates
    obs2.update(yk, uk);

    % Convert dynamic variables to vdata struct
    vars = get_obs_vars(obs2);
    vars_double = {vars.xkp1_est, vars.Pkp1, vars.p_seq_g_Yk, ...
        vars.xkp1_est_f, vars.Pkp1_f};
    vdata = make_data_vectors(vars_double);
    vdata_int16 = make_data_vectors(struct2cell(vars.int16), 'int16');

    % Save observer estimates
    Xk_est{2}(i+1,:) = obs2.xk_est';
    Yk_est{2}(i+1,:) = obs2.yk_est';

end

% Check observer estimates are identical
assert(isequaln(Xk_est{1}, Xk_est{2}))
assert(isequaln(Yk_est{1}, Yk_est{2}))