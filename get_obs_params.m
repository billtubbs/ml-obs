function T = get_obs_params(obs)
% T = get_obs_params(obs) returns a struct containing
% selected parameters of the observer. Which params
% are selected depends on the observer label.

    if strcmp(obs.label, 'none')  % no observer

        T = table();

    elseif startsWith(obs.label, 'KFSS')  % Steady-state Kalman filters

        % Params to return
        params.Q = obs.Q;
        params.R = obs.R;
        T = objects2tablerow(containers.Map({obs.label}, {params}));

    elseif startsWith(obs.label, 'KF')  % Standard Kalman filters

        % Params to return
        params.P0 = obs.P0;
        params.Q = obs.Q;
        params.R = obs.R;
        T = objects2tablerow(containers.Map({obs.label}, {params}));

    elseif startsWith(obs.label, 'LB')  % Luenberger observers

        % Params to return
        params.poles = obs.poles;
        params.K = obs.K;
        T = objects2tablerow(containers.Map({obs.label}, {params}));

    elseif startsWith(obs.label, 'SKF')  % Scheduled Kalman filters

        % Params to return
        params.P0 = obs.P0;
        params.R = obs.R;
        params.Q0 = obs.Q0;
        params.sigma_wp = obs.sigma_wp;
        T = objects2tablerow(containers.Map({obs.label}, {params}));

    elseif startsWith(obs.label, 'RODD')  % RODD MKF observer

        % Params to return
        params.P0 = obs.P0;
        params.Q0 = obs.Q0;
        params.R = obs.R;
        params.epsilon = obs.epsilon;
        params.sigma_wp = obs.sigma_wp;
        params.f = obs.f;
        params.m = obs.m;
        params.d = obs.d;
        params.n_filt = obs.n_filt;
        params.beta = obs.beta;
        T = objects2tablerow(containers.Map({obs.label}, {params}));

    elseif startsWith(obs.label, 'MKF')  % general multi-model Kalman filter

        % Params to return
        %params.P0 = obs.P0;  % Don't store these - too many
        params.Q = obs.Q;
        params.R = obs.R;
        params.f = obs.f;
        params.d = obs.d;
        params.n_filt = obs.n_filt;
        params.beta = obs.beta;
        T = objects2tablerow(containers.Map({obs.label}, {params}));

    elseif startsWith(obs.label, 'AFMM')  % adaptive multi-model Kalman filter

        % Params to return
        params.P0 = obs.P0;
        params.Q0 = obs.Q0;
        params.R = obs.R;
        params.epsilon = obs.epsilon;
        params.sigma_wp = obs.sigma_wp;
        params.n_filt = obs.n_filt;
        params.f = obs.f;
        params.n_min = obs.n_min;
        T = objects2tablerow(containers.Map({obs.label}, {params}));

    elseif startsWith(obs.label, 'EKF')  % Extended Kalman filters

        % Params to return
        params.P0 = obs.P0;
        params.Q = obs.Q;
        params.R = obs.R;
        T = objects2tablerow(containers.Map({obs.label}, {params}));

    elseif startsWith(obs.label, 'MEKF')  % Extended Kalman filters

        % Params to return
        params.P0 = obs.P0;
        params.Q0 = obs.Q0;
        params.R = obs.R;
        params.epsilon = obs.epsilon;
        params.sigma_wp = obs.sigma_wp;
        params.f = obs.f;
        params.m = obs.m;
        params.d = obs.d;
        params.n_filt = obs.n_filt;
        params.beta = obs.beta;
        T = objects2tablerow(containers.Map({obs.label}, {params}));

    else
        error('Value error: observer type not recognized')
    end

end