function T = get_obs_params(obs)
% T = get_obs_params(obs) returns a struct containing
% selected parameters of the observer. Which params
% are selected depends on the observer label.

    switch obs.label
        case 'none'  % no observer
            T = table();

        case {'KF1', 'KF2', 'KF3'}  % Standard Kalman filters

            % Params to return
            params.P0 = obs.P0;
            params.Q = obs.Q;
            params.R = obs.R;
            T = objects2tablerow(containers.Map({obs.label}, {params}));

        case {'LB1', 'LB2'}  % Luenberger observers
            
            % Params to return
            params.poles = obs.poles;
            params.K = obs.K;
            T = objects2tablerow(containers.Map({obs.label}, {params}));
        
        case {'KFSS1', 'KFSS2'}  % Steady-state Kalman filters

            % Params to return
            params.Q = obs.Q;
            params.R = obs.R;
            T = objects2tablerow(containers.Map({obs.label}, {params}));

        case {'SKF', 'SKF1', 'SKF2'}  % Scheduled Kalman filters

            % Params to return
            params.P0 = obs.P0;
            params.R = obs.R;
            params.Q0 = obs.Q0;
            params.sigma_wp = obs.sigma_wp;
            T = objects2tablerow(containers.Map({obs.label}, {params}));

        case {'MKF', 'MKF1', 'MKF2', 'MKF1m', 'MKF2m'}  % RODD MKF observer

            % Params to return
            params.P0 = obs.P0;
            params.Q0 = obs.Q0;
            params.R = obs.R;
            params.epsilon = obs.epsilon;
            params.f = obs.f;
            params.m = obs.m;
            params.d = obs.d;
            params.n_filt = obs.n_filt;
            params.beta = obs.beta;
            T = objects2tablerow(containers.Map({obs.label}, {params}));

        case {'MKF3'}  % general multi-model Kalman filter

            % Params to return
            %params.P0 = obs.P0;  % Don't store these - too many
            params.Q = obs.Q;
            params.R = obs.R;
            params.f = obs.f;
            params.d = obs.d;
            params.n_filt = obs.n_filt;
            params.beta = obs.beta;
            T = objects2tablerow(containers.Map({obs.label}, {params}));
            
        otherwise
            error('Value error: observer type not recognized')
end