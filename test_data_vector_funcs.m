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


%% Test with Kalman filter struct

% Load system and disturbance model from file
sys_rodin_step

% RODD random variable parameters
epsilon = 0.01;
sigma_w = [0.01; 1];

% Sequenece length
nT = 20;

% Simulate system
X0 = zeros(n,1);
t = Ts*(0:nT)';
Wp = sample_random_shocks(nT+1, epsilon, sigma_w(2), sigma_w(1));
U = zeros(nT+1,1);
U(t>=5) = 1;
Y = lsim(Gpss,[U Wp],t,X0);
V = sigma_M*randn(nT+1, 1);
Ym = Y + V;  % measurement

% Load observers from file
obs_rodin_step

% Make copies of KF1
obs1 = KF1;
obs2 = KF1;

% Convert dynamic variables to vdata struct
vars = get_obs_vars(obs2);
vdata = make_data_vectors(struct2cell(vars));

% Simulate observer
Xk_est = {nan(nT+1,n), nan(nT+1,n)};
Yk_est = {nan(nT+1,ny), nan(nT+1,ny)};
for i = 1:nT
    uk = U(i,:)';
    yk = Ym(i,:)';
    obs1 = update_KF(obs1, uk, yk);
    Xk_est{1}(i+1,:) = obs1.xkp1_est';
    Yk_est{1}(i+1,:) = obs1.ykp1_est';

    % Unpack vdata struct and reconstruct observer from KF1
    vars = unpack_data_vectors(vdata);
    assert(isequal(vars{1}, obs2.xkp1_est))
    assert(isequal(vars{2}, obs2.ykp1_est))
    assert(isequal(vars{3}, obs2.P))
    obs2 = KF1;  % makes a new copy
    obs2.xkp1_est(:) = vars{1};
    obs2.ykp1_est(:) = vars{2};
    obs2.P(:) = vars{3};

    % Observer updates
    obs2 = update_KF(obs2, uk, yk);

    % Convert dynamic variables to vdata struct
    vars = get_obs_vars(obs2);
    vdata = make_data_vectors(struct2cell(vars));

    Xk_est{2}(i+1,:) = obs2.xkp1_est';
    Yk_est{2}(i+1,:) = obs2.ykp1_est';
end

% Check observer estimates are identical
assert(isequaln(Xk_est{1}, Xk_est{2}))
assert(isequaln(Yk_est{1}, Yk_est{2}))


%% Test with AFMM MKF observer struct

% Load system and disturbance model from file
sys_rodin_step

% RODD random variable parameters
epsilon = 0.01;
sigma_w = [0.01; 1];

% Sequenece length
nT = 20;

% Simulate system
X0 = zeros(n,1);
t = Ts*(0:nT)';
Wp = sample_random_shocks(nT+1, epsilon, sigma_w(2), sigma_w(1));
U = zeros(nT+1,1);
U(t>=5) = 1;
Y = lsim(Gpss,[U Wp],t,X0);
V = sigma_M*randn(nT+1, 1);
Ym = Y + V;  % measurement

% Load observers from file
obs_rodin_step

% Make copies of AFMM1
obs1 = AFMM1;
obs2 = AFMM1;

% Convert dynamic variables to vdata struct
vars = get_obs_vars(obs2);
vars_double = {vars.xkp1_est, vars.ykp1_est, vars.p_seq_g_Yk, ...
    vars.gamma_k, vars.xkp1_est_f, vars.ykp1_est_f, vars.P_f};
vdata = make_data_vectors(vars_double);
vdata_int16 = make_data_vectors(struct2cell(vars.int16), 'int16');

% Simulate observer
Xk_est = {nan(nT+1,n), nan(nT+1,n)};
Yk_est = {nan(nT+1,ny), nan(nT+1,ny)};
for i = 1:nT
    uk = U(i,:)';
    yk = Ym(i,:)';
    obs1 = update_AFMM(obs1, uk, yk);
    Xk_est{1}(i+1,:) = obs1.xkp1_est';
    Yk_est{1}(i+1,:) = obs1.ykp1_est';

    % Unpack vdata struct and copy values back to observer
    vars_double = unpack_data_vectors(vdata);
    vars_int16 = unpack_data_vectors(vdata_int16);
    vars = struct();
    vars.xkp1_est = vars_double{1};
    vars.ykp1_est = vars_double{2};
    vars.p_seq_g_Yk = vars_double{3};
    vars.gamma_k = vars_double{4};
    vars.xkp1_est_f = vars_double{5};
    vars.ykp1_est_f = vars_double{6};
    vars.P_f = vars_double{7};
    vars.int16.i = vars_int16{1};
    vars.int16.i_next = vars_int16{2};
    vars.int16.f_main = vars_int16{3};
    vars.int16.f_hold = vars_int16{4};
    vars.int16.f_unused = vars_int16{5};
    vars.int16.seq = vars_int16{6};
    obs2_restored = set_obs_vars(AFMM1, vars);  % makes a new copy
    assert(isequal(obs2_restored.xkp1_est, obs2.xkp1_est))
    assert(isequal(obs2_restored.ykp1_est, obs2.ykp1_est))
    assert(isequal(obs2_restored.filters{1}.xkp1_est, obs2.filters{1}.xkp1_est))
    assert(isequal(obs2_restored.filters{1}.ykp1_est, obs2.filters{1}.ykp1_est))
    assert(isequal(obs2_restored.filters{1}.P, obs2.filters{1}.P))
    assert(isequal(obs2_restored.i, obs2.i))
    assert(isequal(obs2_restored.i_next, obs2.i_next))
    assert(isequal(obs2_restored.f_main, obs2.f_main))
    assert(isequal(obs2_restored.f_hold, obs2.f_hold))
    assert(isequal(obs2_restored.f_unused, obs2.f_unused))

    % Observer updates
    obs2 = update_AFMM(obs2, uk, yk);

    % Convert dynamic variables to vdata struct
    vars = get_obs_vars(obs2);
    vars_double = {vars.xkp1_est, vars.ykp1_est, vars.p_seq_g_Yk, ...
        vars.gamma_k, vars.xkp1_est_f, vars.ykp1_est_f, vars.P_f};
    vdata = make_data_vectors(vars_double);
    vdata_int16 = make_data_vectors(struct2cell(vars.int16), 'int16');

    Xk_est{2}(i+1,:) = obs2.xkp1_est';
    Yk_est{2}(i+1,:) = obs2.ykp1_est';
end

% Check observer estimates are identical
assert(isequaln(Xk_est{1}, Xk_est{2}))
assert(isequaln(Yk_est{1}, Yk_est{2}))


%% Test with RODD MKF observer struct

% Load system and disturbance model from file
sys_rodin_step

% RODD random variable parameters
epsilon = 0.01;
sigma_w = [0.01; 1];

% Sequenece length
nT = 20;

% Simulate system
X0 = zeros(n,1);
t = Ts*(0:nT)';
Wp = sample_random_shocks(nT+1, epsilon, sigma_w(2), sigma_w(1));
U = zeros(nT+1,1);
U(t>=5) = 1;
[Y,T,X] = lsim(Gpss,[U Wp],t,X0);
V = sigma_M*randn(nT+1, 1);
Ym = Y + V;  % measurement

% Load observers from file
obs_rodin_step

% Make copies of AFMM1
obs1 = MKF1;
obs2 = MKF1;

% Convert dynamic variables to vdata struct
vars = get_obs_vars(obs2);
vars_double = {vars.xkp1_est, vars.ykp1_est, vars.p_seq_g_Yk, ...
    vars.gamma_k, vars.xkp1_est_f, vars.ykp1_est_f, vars.P_f};
vdata = make_data_vectors(vars_double);
vdata_int16 = make_data_vectors(struct2cell(vars.int16), 'int16');

% Simulate observer
Xk_est = {nan(nT+1,n), nan(nT+1,n)};
Yk_est = {nan(nT+1,ny), nan(nT+1,ny)};
for i = 1:nT
    uk = U(i,:)';
    yk = Ym(i,:)';
    obs1 = update_MKF(obs1, uk, yk);
    Xk_est{1}(i+1,:) = obs1.xkp1_est';
    Yk_est{1}(i+1,:) = obs1.ykp1_est';

    % Unpack vdata struct and copy values back to observer
    vars_double = unpack_data_vectors(vdata);
    vars_int16 = unpack_data_vectors(vdata_int16);
    vars = struct();
    vars.xkp1_est = vars_double{1};
    vars.ykp1_est = vars_double{2};
    vars.p_seq_g_Yk = vars_double{3};
    vars.gamma_k = vars_double{4};
    vars.xkp1_est_f = vars_double{5};
    vars.ykp1_est_f = vars_double{6};
    vars.P_f = vars_double{7};
    vars.int16.i = vars_int16{1};
    vars.int16.i_next = vars_int16{2};
    obs2_restored = set_obs_vars(MKF1, vars);  % makes a new copy
    assert(isequal(obs2_restored.xkp1_est, obs2.xkp1_est))
    assert(isequal(obs2_restored.ykp1_est, obs2.ykp1_est))
    assert(isequal(obs2_restored.filters{1}.xkp1_est, obs2.filters{1}.xkp1_est))
    assert(isequal(obs2_restored.filters{1}.ykp1_est, obs2.filters{1}.ykp1_est))
    assert(isequal(obs2_restored.filters{1}.P, obs2.filters{1}.P))
    assert(isequal(obs2_restored.i, obs2.i))
    assert(isequal(obs2_restored.i_next, obs2.i_next))

    % Observer updates
    obs2 = update_MKF(obs2, uk, yk);

    % Convert dynamic variables to vdata struct
    vars = get_obs_vars(obs2);
    vars_double = {vars.xkp1_est, vars.ykp1_est, vars.p_seq_g_Yk, ...
        vars.gamma_k, vars.xkp1_est_f, vars.ykp1_est_f, vars.P_f};
    vdata = make_data_vectors(vars_double);
    vdata_int16 = make_data_vectors(struct2cell(vars.int16), 'int16');

    Xk_est{2}(i+1,:) = obs2.xkp1_est';
    Yk_est{2}(i+1,:) = obs2.ykp1_est';
end

% Check observer estimates are identical
assert(isequaln(Xk_est{1}, Xk_est{2}))
assert(isequaln(Yk_est{1}, Yk_est{2}))