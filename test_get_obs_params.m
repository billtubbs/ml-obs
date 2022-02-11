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

params = get_obs_params(SKF);
assert(isequal(fieldnames(params), {'P0', 'R', 'Q0', 'sigma_wp'}'))

params = get_obs_params(MKF1);
assert(isequal(fieldnames(params), {'P0', 'Q0', 'R', 'epsilon', 'sigma_wp', 'f', 'm', 'd', 'n_filt', 'beta'}'))
assert(isequal(params.P0, MKF1.P0))
assert(isequal(params.Q0, MKF1.Q0))
assert(isequal(params.R, MKF1.R))
assert(isequal(params.epsilon, MKF1.epsilon))
assert(isequal(params.sigma_wp, MKF1.sigma_wp))
assert(isequal(params.f, MKF1.f))
assert(isequal(params.m, MKF1.m))
assert(isequal(params.d, MKF1.d))
assert(isequal(params.n_filt, MKF1.n_filt))
assert(isequal(params.beta, MKF1.beta))

params = get_obs_params(AFMM1);
assert(isequal(fieldnames(params), {'P0', 'Q0', 'R', 'epsilon', ...
    'sigma_wp', 'n_filt', 'd', 'f', 'n_min'}'))
assert(isequal(params.P0, AFMM1.P0))
assert(isequal(params.Q0, AFMM1.Q0))
assert(isequal(params.R, AFMM1.R))
assert(isequal(params.epsilon, AFMM1.epsilon))
assert(isequal(params.sigma_wp, AFMM1.sigma_wp))
assert(isequal(params.n_filt, AFMM1.n_filt))
assert(isequal(params.d, AFMM1.d))
assert(isequal(params.f, AFMM1.f))
assert(isequal(params.n_min, AFMM1.n_min))
