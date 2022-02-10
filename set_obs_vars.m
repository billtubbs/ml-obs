function obs = set_obs_vars(obs, vars)
% obs = set_obs_vars(obs, vars) sets all the variables (i.e.
% time-varying properties) of the observer from the values
% in vars (a struct).
%
    if strcmp(obs.label, 'none')  % no observer

        % No variables

    elseif startsWith(obs.label, 'KFSS')  % Steady-state Kalman filters

        % Set variables
        obs.xkp1_est = vars.xkp1_est;
        obs.ykp1_est = vars.ykp1_est;

    elseif startsWith(obs.label, 'KF')  % Standard Kalman filters

        % Set variables
        obs.xkp1_est = vars.xkp1_est;
        obs.ykp1_est = vars.ykp1_est;
        obs.P = vars.P;

    elseif startsWith(obs.label, 'LB')  % Luenberger observers

        % Set variables
        obs.xkp1_est = vars.xkp1_est;
        obs.ykp1_est = vars.ykp1_est;

    elseif startsWith(obs.label, 'SKF')  % Scheduled Kalman filters

        % Set variables
        obs.xkp1_est = vars.xkp1_est;
        obs.ykp1_est = vars.ykp1_est;
        obs.P = vars.P;

    elseif startsWith(obs.label, 'MMKF')  % general multi-model Kalman filter

        % Set variables
        obs.xkp1_est = vars.xkp1_est;
        obs.ykp1_est = vars.ykp1_est;
        for f = 1:obs.n_filt
           obs.filters{f}.xkp1_est = vars.xkp1_est_f{f};
           obs.filters{f}.P = vars.P_f{f};
        end

    elseif startsWith(obs.label, 'MKF')  % RODD MKF observer

        % Set double variables
        obs.xkp1_est = vars.xkp1_est;
        obs.ykp1_est = vars.ykp1_est;
        obs.p_seq_g_Yk = vars.p_seq_g_Yk;
        for f = 1:obs.n_filt
           obs.filters{f}.xkp1_est = vars.xkp1_est_f{f};
           obs.filters{f}.ykp1_est = vars.ykp1_est_f{f};
           obs.filters{f}.P = vars.P_f{f};
        end
        % Set integer variables
        obs.i = vars.int16.i;
        obs.i_next = vars.int16.i_next;

        % TODO: Are any of these others dynamic?
        % vars.p_yk_g_seq_Ykm1 = obs.p_yk_g_seq_Ykm1;
        % vars.p_gammak_g_Ykm1 = obs.p_gammak_g_Ykm1;
        % vars.p_gamma_k = obs.p_gamma_k;
        % vars.p_seq_g_Ykm1 = obs.p_seq_g_Ykm1;

    elseif startsWith(obs.label, 'AFMM')  % adaptive multi-model Kalman filter

        % Set double variables
        obs.xkp1_est = vars.xkp1_est;
        obs.ykp1_est = vars.ykp1_est;
        obs.p_seq_g_Yk = vars.p_seq_g_Yk;
        for f = 1:obs.n_filt
           obs.filters{f}.xkp1_est = vars.xkp1_est_f{f};
           obs.filters{f}.ykp1_est = vars.ykp1_est_f{f};
           obs.filters{f}.P = vars.P_f{f};
        end
        % Set integer variables
        obs.i = vars.int16.i;
        obs.i_next = vars.int16.i_next;
        obs.f_main = vars.int16.f_main;
        obs.f_hold = vars.int16.f_hold;
        obs.f_unused = vars.int16.f_unused;
        obs.seq = vars.int16.seq;

        % TODO: Are any of these others dynamic?
        % vars.p_yk_g_seq_Ykm1 = obs.p_yk_g_seq_Ykm1;
        % vars.p_gammak_g_Ykm1 = obs.p_gammak_g_Ykm1;
        % vars.p_gamma_k = obs.p_gamma_k;
        % vars.p_seq_g_Ykm1 = obs.p_seq_g_Ykm1;

    elseif startsWith(obs.label, 'EKF')  % Extended Kalman filters

        % Set variables
        % TODO: Add dynamic vars

    elseif startsWith(obs.label, 'MEKF')  % Extended Kalman filters

        % Set variables
        % TODO: Add dynamic vars

    else
        error('Value error: observer type not recognized')
    end

end