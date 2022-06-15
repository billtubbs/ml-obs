function params = get_obs_params(obs)
% params = get_obs_params(obs) returns a struct containing
% selected parameters of the observer. Which params
% are selected depends on the observer type.
%
    switch(obs.type)
        
        case 'none'  % no observer

            params = struct();  % No parameters

        case 'KFSS'  % Steady-state Kalman filters

            % Params to return
            params.Q = obs.Q;
            params.R = obs.R;

        case 'KF'  % Standard Kalman filters

            % Params to return
            params.P0 = obs.P0;
            params.Q = obs.Q;
            params.R = obs.R;

        case 'LB'  % Luenberger observers

            % Params to return
            params.poles = obs.poles;
            params.K = obs.K;

        case 'SKF'  % Scheduled Kalman filters

            % Params to return
            params.P0 = obs.P0;  % Don't store these - too many
            params.Q = obs.Q;
            params.R = obs.R;
            params.f = obs.f;
            params.d = obs.d;
            params.n_filt = obs.n_filt;

        case 'MKF'  % general multi-model Kalman filter

            % Params to return
            params.P0 = obs.P0;  % Don't store these - too many
            params.Q = obs.Q;
            params.R = obs.R;
            params.f = obs.f;
            params.d = obs.d;
            params.n_filt = obs.n_filt;
            params.beta = obs.beta;

        case 'MKF_SF'  % RODD MKF observer

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

        case 'MKF_SP'  % adaptive multi-model Kalman filter

            % Params to return
            params.P0 = obs.P0;
            params.Q0 = obs.Q0;
            params.R = obs.R;
            params.epsilon = obs.epsilon;
            params.sigma_wp = obs.sigma_wp;
            params.n_filt = obs.n_filt;
            params.d = obs.d;
            params.f = obs.f;
            params.n_min = obs.n_min;

        case 'EKF'  % Extended Kalman filters

            % Params to return
            params.P0 = obs.P0;
            params.Q = obs.Q;
            params.R = obs.R;

        case 'MEKF'  % Extended Kalman filters

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

        otherwise
            error('Value error: observer type not recognized')
    end

end