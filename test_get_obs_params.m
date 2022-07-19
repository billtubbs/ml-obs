% test get_obs_params.m

clear all

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step


%% Test all fieldnames returned

params = get_obs_params(KFSS1);
assert(isequal(fieldnames(params), {'Q', 'R'}'))
assert(isequal(params.Q, KFSS1.Q))
assert(isequal(params.R, KFSS1.R))

params = get_obs_params(KF1);
assert(isequal(fieldnames(params), {'P0', 'Q', 'R'}'))
assert(isequal(params.P0, KF1.P0))

params = get_obs_params(LB1);
assert(isequal(fieldnames(params), {'poles', 'K'}'))
assert(isequal(params.poles, LB1.poles))

% Define a scheduled Kalman filter
nT = 100;
t = Ts*(0:nT)';
A2 = repmat({A}, 1, 2);
Bu2 = repmat({Bu}, 1, 2);
C2 = repmat({C}, 1, 2);
Du2 = repmat({Du}, 1, 2);
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
P0_init = repmat({P0}, 1, 2);
Q2 = {diag([Q0(1,1) sigma_wp(1,1)^2]), ...
      diag([Q0(1,1) sigma_wp(1,2)^2])};
R2 = {sigma_M^2, sigma_M^2};
seq = {zeros(1, nT+1)};
seq{1}(t == 10) = 1;
SKF = MKFObserverSched(A2,Bu2,C2,Ts,P0,Q2,R2,seq{1},"SKF");

params = get_obs_params(SKF);
assert(isequal(fieldnames(params), {'P0', 'Q', 'R', 'f', 'd', 'n_filt'}'))

params = get_obs_params(MKF3);
assert(isequal(fieldnames(params), {'P0', 'Q', 'R', 'f', 'n_filt'}'))

params = get_obs_params(MKF_SF1);
assert(isequal(fieldnames(params), {'P0', 'Q0', 'R', 'epsilon', ...
    'sigma_wp', 'f', 'm', 'd', 'n_filt', 'beta'}'))
assert(isequal(params.P0, MKF_SF1.P0))
assert(isequal(params.Q0, MKF_SF1.Q0))
assert(isequal(params.R, MKF_SF1.R))
assert(isequal(params.epsilon, MKF_SF1.epsilon))
assert(isequal(params.sigma_wp, MKF_SF1.sigma_wp))
assert(isequal(params.f, MKF_SF1.f))
assert(isequal(params.m, MKF_SF1.m))
assert(isequal(params.d, MKF_SF1.d))
assert(isequal(params.n_filt, MKF_SF1.n_filt))
assert(isequal(params.beta, MKF_SF1.beta))

params = get_obs_params(MKF_SF95);
assert(isequal(fieldnames(params), {'P0', 'Q0', 'R', 'epsilon', ...
    'sigma_wp', 'f', 'm', 'd', 'n_filt', 'beta'}'))
assert(isequal(params.P0, MKF_SF95.P0))
assert(isequal(params.Q0, MKF_SF95.Q0))
assert(isequal(params.R, MKF_SF95.R))
assert(isequal(params.epsilon, MKF_SF95.epsilon))
assert(isequal(params.sigma_wp, MKF_SF95.sigma_wp))
assert(isequal(params.f, MKF_SF95.f))
assert(isequal(params.m, MKF_SF95.m))
assert(isequal(params.d, MKF_SF95.d))
assert(isequal(params.n_filt, MKF_SF95.n_filt))
assert(isequal(params.beta, MKF_SF95.beta))

params = get_obs_params(MKF_SP1);
assert(isequal(fieldnames(params), {'P0', 'Q0', 'R', 'epsilon', ...
    'sigma_wp', 'n_filt', 'f', 'n_min'}'))
assert(isequal(params.P0, MKF_SP1.P0))
assert(isequal(params.Q0, MKF_SP1.Q0))
assert(isequal(params.R, MKF_SP1.R))
assert(isequal(params.epsilon, MKF_SP1.epsilon))
assert(isequal(params.sigma_wp, MKF_SP1.sigma_wp))
assert(isequal(params.n_filt, MKF_SP1.n_filt))
assert(isequal(params.f, MKF_SP1.f))
assert(isequal(params.n_min, MKF_SP1.n_min))
