% test get_obs_params.m

clear all

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step


%% Test all fieldnames returned

params = get_obs_params(KFSS1);
assert(isequal(fieldnames(params), {'Q', 'R'}'))

params = get_obs_params(KF1);
assert(isequal(fieldnames(params), {'P0', 'Q', 'R'}'))

params = get_obs_params(LB1);
assert(isequal(fieldnames(params), {'poles', 'K'}'))

params = get_obs_params(SKF);
assert(isequal(fieldnames(params), {'P0', 'R', 'Q0', 'sigma_wp'}'))

params = get_obs_params(MKF1);
assert(isequal(fieldnames(params), {'Q', 'R', 'f', 'd', 'n_filt', 'beta'}'))

% TODO: Rename the MKF1, MKF2 observers
%params = get_obs_params(RODD);
%assert(isequal(fieldnames(params), {'Q', 'R', 'f', 'd', 'n_filt', 'beta'}'))

params = get_obs_params(AFMM1);
assert(isequal(fieldnames(params), {'P0', 'Q0', 'R', 'epsilon', ...
    'sigma_wp', 'n_filt', 'f', 'n_min'}'))
