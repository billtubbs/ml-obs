function obs = set_obs_vars(obs, vars)
% obs = set_obs_vars(obs, vars) sets all the variables (i.e.
% time-varying properties) of the observer from the values
% in vars (a struct).
%

% TODO: Write a test script for set_obs_vars and get_obs_vars
    switch(obs.type)

        case {"LB", "KFPSS", "KFFSS"}  % steady-state filters

            % Set variables
            obs.xkp1_est = vars.xkp1_est;
            if ismember(obs.type, ["LB" "KFPSS"])
                obs.ykp1_est = vars.ykp1_est;
            end

        case "KFP"  % Kalman filter - prediction form

            % Set variables
            obs.xkp1_est = vars.xkp1_est;
            obs.ykp1_est = vars.ykp1_est;
            obs.Pkp1 = vars.Pkp1;

        case "KFF"  % Kalman filter - filtering form

            % Set variables
            obs.xkp1_est = vars.xkp1_est;
            obs.Pkp1 = vars.Pkp1;

        case "SKF"  % Switching Kalman filter

            % Set variables
            obs.xkp1_est = vars.xkp1_est;
            obs.Pkp1 = vars.Pkp1;
            obs.rk = vars.rk;

        case "SKF_S"  % Switching Kalman filter with sequence

            % Set variables
            obs.xkp1_est = vars.xkp1_est;
            obs.Pkp1 = vars.Pkp1;

            % Set integer variables
            if startsWith(obs.type, "SKF_S")
                obs.i = vars.int16.i;
                obs.i_next = vars.int16.i_next;
            end

        case {"MKF", "MKF_SF", "MKF_SF_RODD95"}  % Multi-model Kalman filters

            % Set double variables
            obs.xkp1_est = vars.xkp1_est;
            obs.Pkp1 = vars.Pkp1;
            obs.p_seq_g_Yk = vars.p_seq_g_Yk;
            obs.filters.Xkp1_est = vars.xkp1_est_f;
            obs.filters.Pkp1 = vars.Pkp1_f;

            % Set integer variables
            obs.rk = vars.int16.rk;
            if strcmp(obs.type, "MKF_S")
                obs.i = vars.int16.i;
                obs.i_next = vars.int16.i_next;
                obs.seq = vars.int16.seq;
            end

        case {"MKF_DI", "MKF_SF_RODD"}

            % Set double variables
            obs.xkp1_est = vars.xkp1_est;
            obs.p_seq_g_Yk = vars.p_seq_g_Yk;
            obs.filters.Xkp1_est = vars.xkp1_est_f;
            obs.filters.Pkp1 = vars.Pkp1_f;

            % Set integer variables
            obs.rk = vars.int16.rk;
            obs.i = vars.int16.i;
            obs.i_next = vars.int16.i_next;
            obs.i2 = vars.int16.i2;
            obs.i2_next = vars.int16.i2_next;

        case "MKF_SP"  % MKF observer with sequence pruning

            % Set double variables
            obs.xkp1_est = vars.xkp1_est;
            obs.p_seq_g_Yk = vars.p_seq_g_Yk;
            obs.filters.Xkp1_est = vars.xkp1_est_f;
            obs.filters.Pkp1 = vars.Pkp1_f;

            % Set integer variables
            obs.rk = vars.int16.rk;
            obs.f_main = vars.int16.f_main;
            obs.f_hold = vars.int16.f_hold;

        case "EKF"  % Extended Kalman filters

            % Set variables
            % TODO: Add dynamic vars

        case "MEKF"  % Extended Kalman filters

            % Set variables
            % TODO: Add dynamic vars

        otherwise
            error("Value error: observer type not recognized")
    end

end