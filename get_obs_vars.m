function vars = get_obs_vars(obs)
% vars = get_obs_vars(obs) returns a struct containing
% all the variables (i.e. time-varying properties) of the
% observer. The variables returned depends on the observer type.
%
    if strcmp(obs.label, 'none')  % no observer

        vars = struct();

    elseif startsWith(obs.label, 'KFSS')  % Steady-state Kalman filters

        % Vars to return
        vars.xkp1_est = obs.xkp1_est;
        vars.ykp1_est = obs.ykp1_est;

    elseif startsWith(obs.label, 'KF')  % Standard Kalman filters

        % Vars to return
        vars.xkp1_est = obs.xkp1_est;
        vars.ykp1_est = obs.ykp1_est;
        vars.P = obs.P;

    elseif startsWith(obs.label, 'LB')  % Luenberger observers

        % Vars to return
        vars.xkp1_est = obs.xkp1_est;
        vars.ykp1_est = obs.ykp1_est;

    elseif startsWith(obs.label, 'SKF')  % Scheduled Kalman filters

        % Vars to return (same as Kalman filter)
        vars.xkp1_est = obs.xkp1_est;
        vars.ykp1_est = obs.ykp1_est;
        vars.P = obs.P;

    elseif startsWith(obs.label, 'MKF')  % general multi-model Kalman filter

        % Vars to return
        vars.xkp1_est = obs.xkp1_est;
        vars.ykp1_est = obs.ykp1_est;
        vars.P_f = cell(1, obs.n_filt);
        vars.xkp1_est_f = cell(1, obs.n_filt);
        for f = 1:obs.n_filt
           vars.xkp1_est_f{f} = obs.filters{f}.xkp1_est;
           vars.P_f{f} = obs.filters{f}.P;
        end

    elseif startsWith(obs.label, 'RODD')  % RODD MKF observer

        % Vars to return
        vars.xkp1_est = obs.xkp1_est;
        vars.ykp1_est = obs.ykp1_est;
        vars.xkp1_est = obs.xkp1_est;
        vars.ykp1_est = obs.ykp1_est;
        vars.p_seq_g_Yk = obs.p_seq_g_Yk;
        vars.P_f = cell(1, obs.n_filt);
        vars.xkp1_est_f = cell(1, obs.n_filt);
        for f = 1:obs.n_filt
           vars.xkp1_est_f{f} = obs.filters{f}.xkp1_est;
           vars.P_f{f} = obs.filters{f}.P;
        end
        % Integer variables
        vars.int16.i = obs.i;
        vars.int16.i_next = obs.i_next;
        vars.int16.f_main = obs.f_main;
        vars.int16.f_hold = obs.f_hold;
        vars.int16.f_unused = obs.f_unused;

        % TODO: Are any of these others dynamic?
        % vars.p_yk_g_seq_Ykm1 = obs.p_yk_g_seq_Ykm1;
        % vars.p_gammak_g_Ykm1 = obs.p_gammak_g_Ykm1;
        % vars.p_gamma_k = obs.p_gamma_k;
        % vars.p_seq_g_Ykm1 = obs.p_seq_g_Ykm1;

    elseif startsWith(obs.label, 'AFMM')  % adaptive multi-model Kalman filter

        % Vars to return
        vars.xkp1_est = obs.xkp1_est;
        vars.ykp1_est = obs.ykp1_est;
        vars.xkp1_est = obs.xkp1_est;
        vars.ykp1_est = obs.ykp1_est;
        vars.p_seq_g_Yk = obs.p_seq_g_Yk;
        vars.P_f = cell(1, obs.n_filt);
        vars.xkp1_est_f = cell(1, obs.n_filt);
        for f = 1:obs.n_filt
           vars.xkp1_est_f{f} = obs.filters{f}.xkp1_est;
           vars.P_f{f} = obs.filters{f}.P;
        end
        % Integer variables
        vars.int16.i = obs.i;
        vars.int16.i_next = obs.i_next;
        vars.int16.f_main = obs.f_main;
        vars.int16.f_hold = obs.f_hold;
        vars.int16.f_unused = obs.f_unused;
        vars.int16.seq = obs.seq;

        % TODO: Are any of these others dynamic?
        % vars.p_yk_g_seq_Ykm1 = obs.p_yk_g_seq_Ykm1;
        % vars.p_gammak_g_Ykm1 = obs.p_gammak_g_Ykm1;
        % vars.p_gamma_k = obs.p_gamma_k;
        % vars.p_seq_g_Ykm1 = obs.p_seq_g_Ykm1;

    elseif startsWith(obs.label, 'EKF')  % Extended Kalman filters

        % Vars to return
        % TODO: Add dynamic vars

    elseif startsWith(obs.label, 'MEKF')  % Extended Kalman filters

        % Vars to return
        % TODO: Add dynamic vars

    else
        error('Value error: observer type not recognized')
    end

end