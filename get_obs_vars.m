function vars = get_obs_vars(obs)
% vars = get_obs_vars(obs) returns a struct containing
% all the variables (i.e. time-varying properties) of the
% observer. The variables returned depends on the observer type.
%

% TODO: Write a test script for set_obs_vars and get_obs_vars
    switch(obs.type)

        case {"KFSS", "LB"}  % Steady-state filters

            % Vars to return
            vars.xkp1_est = obs.xkp1_est;

        case {"KFF", "SKF", "SKF_S"}  % Kalman filter and switching KF

            % Vars to return
            vars.xkp1_est = obs.xkp1_est;
            vars.Pkp1 = obs.Pkp1;

        case {"MKF", "MKF_S", "MKF_SF", "MKF_SP"}  % multi-model Kalman filters

            % Vars to return
            vars.xkp1_est = obs.xkp1_est;
            vars.Pkp1 = obs.Pkp1;
            vars.p_seq_g_Yk = obs.p_seq_g_Yk;
            vars.rk = obs.rk;
            vars.xkp1_est_f = obs.filters.Xkp1_est;
            vars.Pkp1_f = obs.filters.Pkp1_est;

            if startsWith(obs.type, "MKF_SP")
                % Additional variables used by adaptive sequence
                % pruning algorithms
                vars.int16.f_main = obs.f_main;
                vars.int16.f_hold = obs.f_hold;
            end

        case {"MKF_SFF", "MKF_SPF"}  % special multi-model Kalman filters

            % Vars to return
            vars.xkp1_est = obs.xkp1_est;
            vars.Pkp1 = obs.Pkp1;
            vars.p_seq_g_Yk = obs.p_seq_g_Yk;
            vars.gamma_k = obs.gamma_k;
            vars.xkp1_est_f = cell(1, obs.n_filt);
            vars.ykp1_est_f = cell(1, obs.n_filt);
            vars.Pkp1_f = cell(1, obs.n_filt);
            for f = 1:obs.n_filt
               vars.xkp1_est_f{f} = obs.filters{f}.xkp1_est;
               vars.ykp1_est_f{f} = obs.filters{f}.ykp1_est;
               vars.Pkp1_f{f} = obs.filters{f}.Pkp1;
            end
            % Integer variables
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

            % Vars to return
            % TODO: Add dynamic vars

        case "MEKF"  % Extended Kalman filters

            % Vars to return
            % TODO: Add dynamic vars

        otherwise
            error("Value error: observer type not recognized")
    end

end