function params = get_obs_params(obs)
% params = get_obs_params(obs) returns a struct containing
% selected parameters of the observer. Which params
% are selected depends on the observer type.
%
    switch(obs.type)

        case "KFSS"  % Steady-state Kalman filters

            % Params to return
            params.A = obs.A;
            params.B = obs.B;
            params.C = obs.C;
            params.Q = obs.Q;
            params.R = obs.R;

        case "KF"  % Standard Kalman filters

            % Params to return
            params.A = obs.A;
            params.B = obs.B;
            params.C = obs.C;
            params.Q = obs.Q;
            params.R = obs.R;
            params.P0 = obs.P0;

        case "KFF"  % New Kalman filters

            % Params to return
            params.model = obs.model;
            params.P0 = obs.P0;

        case "LB"  % Luenberger observers

            % Params to return
            params.A = obs.A;
            params.B = obs.B;
            params.C = obs.C;
            params.poles = obs.poles;
            params.K = obs.K;

        case {"SKF", "SKF_S"}  % Switching Kalman filters

            % Params to return
            params.models = obs.models;
            params.P0 = obs.P0;
            params.nj = obs.nj;

        case "MKF"  % general multi-model Kalman filters

            % Params to return
            params.models = obs.models;
            params.P0 = obs.P0;
            params.T = obs.T;
            params.nh = obs.nh;

        case "MKF_DI"  % multi-model Kalman filter with detection interval

            % Params to return
            params.models = obs.models;
            params.P0 = obs.P0;
            params.T = obs.T;
            params.d = obs.d;
            params.nh = obs.nh;
            params.nf = obs.nf;

        case {"MKF_SF"}  % MKF observer with sequence fusion

            % Params to return
            params.models = obs.models;
            params.P0 = obs.P0;
            params.T = obs.T;
            params.nh = obs.nh;
            params.nf = obs.nf;
        
        case {"MKF_SF_RODD", "MKF_SF_RODD95"}  % MKF observers for RODDs

            % Params to return
            params.model = obs.sys_model;
            params.P0 = obs.P0;
            params.Q0 = obs.Q0;
            params.R = obs.R;
            params.epsilon = obs.epsilon;
            params.sigma_wp = obs.sigma_wp;
            params.f = obs.f;
            params.m = obs.m;
            params.d = obs.d;
            params.nh = obs.nh;
            params.beta = obs.beta;

        case {"MKF_SP", "MKF_SPF"}   % MKF observer with sequence pruning

            % Params to return
            params.model = obs.sys_model;
            params.P0 = obs.P0;
            params.Q0 = obs.Q0;
            params.R = obs.R;
            params.epsilon = obs.epsilon;
            params.sigma_wp = obs.sigma_wp;
            params.nh = obs.nh;
            params.n_min = obs.n_min;

        case "EKF"  % Extended Kalman filters

            % Params to return
            params.P0 = obs.P0;
            params.Q = obs.Q;
            params.R = obs.R;

        case "MEKF"  % Extended Kalman filters

            % Params to return
            params.P0 = obs.P0;
            params.T = obs.T;
            params.nh = obs.nh;

        otherwise
            error("Value error: observer type not recognized")
    end

end