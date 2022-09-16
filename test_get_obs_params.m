% test get_obs_params.m

clear all

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step


%% Test all fieldnames returned

params = get_obs_params(KFSS1);
assert(isequal(fieldnames(params), {'A', 'B', 'C', 'Q', 'R'}'))
assert(isequal(params.A, KFSS1.A))
assert(isequal(params.B, KFSS1.B))
assert(isequal(params.C, KFSS1.C))
assert(isequal(params.Q, KFSS1.Q))
assert(isequal(params.R, KFSS1.R))

params = get_obs_params(KF1);
assert(isequal(fieldnames(params), {'model', 'P0'}'))
assert(isequal(params.model, KF1.model))
assert(isequal(params.P0, KF1.P0))

params = get_obs_params(LB1);
assert(isequal(fieldnames(params), {'A', 'B', 'C', 'poles', 'K'}'))
assert(isequal(params.A, LB1.A))
assert(isequal(params.B, LB1.B))
assert(isequal(params.C, LB1.C))
assert(isequal(params.poles, LB1.poles))
assert(isequal(params.K, LB1.K))

nT = 100;
t = Ts*(0:nT)';

% Define a scheduled Kalman filter
P0 = 1000*eye(n);
Q1 = diag([Q0(1,1) sigma_wp(1,1)^2]);
Q2 = diag([Q0(1,1) sigma_wp(1,2)^2]);
R = sigma_M^2;
models = {model, model};
models{1}.Bu = Bu;
models{2}.Bu = Bu;
models{1}.Q = Q1;
models{2}.Q = Q2;
models{1}.R = R;
models{2}.R = R;
seq = {zeros(1, nT+1)};
seq{1}(t == 10) = 1;
SKF1 = SKFObserverS(models,P0,seq{1},"SKF1");

params = get_obs_params(SKF1);
assert(isequal(fieldnames(params), {'models', 'P0', 'nj'}'))

% TODO: test other types of MKF
% params = get_obs_params(MKF3);
% assert(isequal(fieldnames(params), {'P0', 'Q', 'R', 'f', 'nh'}'))

params = get_obs_params(MKF_SF1);
assert(isequal(fieldnames(params), {'model', 'P0', 'Q0', 'R', 'epsilon', ...
    'sigma_wp', 'f', 'm', 'd', 'nh', 'beta'}'))
assert(isequal(params.model, MKF_SF1.sys_model))
assert(isequal(params.P0, MKF_SF1.P0))
assert(isequal(params.Q0, MKF_SF1.Q0))
assert(isequal(params.R, MKF_SF1.R))
assert(isequal(params.epsilon, MKF_SF1.epsilon))
assert(isequal(params.sigma_wp, MKF_SF1.sigma_wp))
assert(isequal(params.f, MKF_SF1.f))
assert(isequal(params.m, MKF_SF1.m))
assert(isequal(params.d, MKF_SF1.d))
assert(isequal(params.nh, MKF_SF1.nh))
assert(isequal(params.beta, MKF_SF1.beta))

% params = get_obs_params(MKF_SF95);
% assert(isequal(fieldnames(params), {'model', 'P0', 'Q0', 'R', 'epsilon', ...
%     'sigma_wp', 'f', 'm', 'd', 'nh', 'beta'}'))
% assert(isequal(params.P0, MKF_SF95.P0))
% assert(isequal(params.Q0, MKF_SF95.Q0))
% assert(isequal(params.R, MKF_SF95.R))
% assert(isequal(params.epsilon, MKF_SF95.epsilon))
% assert(isequal(params.sigma_wp, MKF_SF95.sigma_wp))
% assert(isequal(params.f, MKF_SF95.f))
% assert(isequal(params.m, MKF_SF95.m))
% assert(isequal(params.d, MKF_SF95.d))
% assert(isequal(params.nh, MKF_SF95.nh))
% assert(isequal(params.beta, MKF_SF95.beta))

params = get_obs_params(MKF_SP1);
assert(isequal(fieldnames(params), {'model', 'P0', 'Q0', 'R', 'epsilon', ...
    'sigma_wp', 'nh', 'n_min'}'))
assert(isequal(params.model, MKF_SP1.sys_model))
assert(isequal(params.P0, MKF_SP1.P0))
assert(isequal(params.Q0, MKF_SP1.Q0))
assert(isequal(params.R, MKF_SP1.R))
assert(isequal(params.epsilon, MKF_SP1.epsilon))
assert(isequal(params.sigma_wp, MKF_SP1.sigma_wp))
assert(isequal(params.nh, MKF_SP1.nh))
assert(isequal(params.n_min, MKF_SP1.n_min))
