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
assert(isequal(fieldnames(params), {'Q', 'R', 'f', 'd', 'n_filt', 'beta'}'))
assert(isequal(params.Q, MKF1.Q))
assert(isequal(params.R, MKF1.R))
assert(isequal(params.f, MKF1.f))
assert(isequal(params.d, MKF1.d))
assert(isequal(params.beta, MKF1.beta))

% TODO: Rename the MKF1, MKF2 observers
%params = get_obs_params(RODD);
%assert(isequal(fieldnames(params), {'Q', 'R', 'f', 'd', 'n_filt', 'beta'}'))

params = get_obs_params(AFMM1);
assert(isequal(fieldnames(params), {'P0', 'Q0', 'R', 'epsilon', ...
    'sigma_wp', 'n_filt', 'd', 'f', 'n_min'}'))
assert(isequal(params.Q0, AFMM1.Q0))
assert(isequal(params.P0, AFMM1.P0))
assert(isequal(params.f, AFMM1.f))
assert(isequal(params.d, AFMM1.d))
assert(isequal(params.n_filt, AFMM1.n_filt))