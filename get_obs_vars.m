function vars = get_obs_vars(obs)
% vars = get_obs_vars(obs) returns a struct containing
% all the variables (i.e. time-varying properties) of the
% observer. The variables returned depends on the observer type.
%

% TODO: Write a test script for set_obs_vars and get_obs_vars
    switch(obs.type)

        case {"KFPSS", "LB"}  % Steady-state filters

            % Get variables
            vars.xkp1_est = obs.xkp1_est;
            vars.ykp1_est = obs.ykp1_est;

        case {"KFFSS"}  % Steady-state filters

            % Get variables
            vars.xkp1_est = obs.xkp1_est;

        case {"KFP"}  % Kalman filter prediction

            % Get variables
            vars.xkp1_est = obs.xkp1_est;
            vars.ykp1_est = obs.ykp1_est;
            vars.Pkp1 = obs.Pkp1;

        case {"KFF", "SKF", "SKF_S"}  % Kalman filter and switching KF

            % Get variables
            vars.xkp1_est = obs.xkp1_est;
            vars.Pkp1 = obs.Pkp1;

        case {"MKF_DI", "MKF_SF_RODD", "MKF_SF_RODD95"}

            % Get double variables
            vars.xkp1_est = obs.xkp1_est;
            vars.Pkp1 = obs.Pkp1;
            vars.p_seq_g_Yk = obs.p_seq_g_Yk;
            vars.xkp1_est_f = obs.filters.Xkp1_est;
            vars.Pkp1_f = obs.filters.Pkp1;

            % Get integer variables
            vars.int16.rk = obs.rk;
            vars.int16.i = obs.i;
            vars.int16.i_next = obs.i_next;
            vars.int16.i2 = obs.i2;
            vars.int16.i2_next = obs.i2_next;
            vars.int16.seq = obs.seq;

        case {"MKF", "MKF_S", "MKF_SF", "MKF_SP"}  % multi-model Kalman filters

            % Get double variables
            vars.xkp1_est = obs.xkp1_est;
            vars.Pkp1 = obs.Pkp1;
            vars.p_seq_g_Yk = obs.p_seq_g_Yk;
            vars.xkp1_est_f = obs.filters.Xkp1_est;
            vars.Pkp1_f = obs.filters.Pkp1;

            % Get integer variables
            vars.int16.rk = obs.rk;
            if startsWith(obs.type, "MKF_SP")
                % Additional variables used by adaptive sequence
                % pruning algorithms
                vars.int16.f_main = obs.f_main;
                vars.int16.f_hold = obs.f_hold;
            end

        case {"MKF_SFF", "MKF_SPF"}  % special multi-model Kalman filters

            % Get double variables
            vars.xkp1_est = obs.xkp1_est;
            vars.Pkp1 = obs.Pkp1;
            vars.p_seq_g_Yk = obs.p_seq_g_Yk;
            vars.xkp1_est_f = nan(obs.n, 1, obs.nh);
            vars.Pkp1_f = nan(obs.n, obs.n, obs.nh);
            for f = 1:obs.nh
               vars.xkp1_est_f(:,:,f) = obs.filters{f}.xkp1_est;
               vars.Pkp1_f(:,:,f) = obs.filters{f}.Pkp1;
            end

            % Get integer variables
            vars.int16.rk = obs.rk;
            vars.int16.i = obs.i;
            vars.int16.i_next = obs.i_next;
            if strcmp(obs.type, "MKF_SPF")
                % Additional variables used by adaptive sequence
                % pruning algorithms
                vars.int16.f_main = obs.f_main;
                vars.int16.f_hold = obs.f_hold;
                vars.int16.seq = obs.seq;
            end

        case "EKF"  % Extended Kalman filters

            % Get double variables
            % TODO: Add dynamic vars

        case "MEKF"  % Extended Kalman filters

            % Get double variables
            % TODO: Add dynamic vars

        otherwise
            error("Value error: observer type not recognized")
    end

end