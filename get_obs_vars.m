function vars = get_obs_vars(obs)
% vars = get_obs_vars(obs) returns a struct containing
% all the variables (i.e. time-varying properties) of the
% observer. The variables returned depends on the observer type.
%

% TODO: Write a test script for set_obs_vars and get_obs_vars
    switch(obs.type)

        case 'none'  % no observer

            vars = struct();  % No variables

        case {'KFSS', 'LB'}  % Steady-state filters

            % Vars to return
            vars.xkp1_est = obs.xkp1_est;
            vars.ykp1_est = obs.ykp1_est;

        case 'KF'  % Standard Kalman filters

            % Vars to return
            vars.xkp1_est = obs.xkp1_est;
            vars.ykp1_est = obs.ykp1_est;
            vars.P = obs.P;

        case 'SKF'  % Scheduled Kalman filters

            % Vars to return (same as Kalman filter)
            vars.xkp1_est = obs.xkp1_est;
            vars.ykp1_est = obs.ykp1_est;
            vars.P = obs.P;

        case 'MKF'  % general multi-model Kalman filter

            % Vars to return
            vars.xkp1_est = obs.xkp1_est;
            vars.ykp1_est = obs.ykp1_est;
            vars.P_f = cell(1, obs.n_filt);
            vars.ykp1_est_f = cell(1, obs.n_filt);
            vars.xkp1_est_f = cell(1, obs.n_filt);
            for f = 1:obs.n_filt
               vars.xkp1_est_f{f} = obs.filters{f}.xkp1_est;
               vars.ykp1_est_f{f} = obs.filters{f}.ykp1_est;
               vars.P_f{f} = obs.filters{f}.P;
            end

        case {"MKF-SP", "MKF-SF"}  % multi-model Kalman filters

            % Vars to return
            vars.xkp1_est = obs.xkp1_est;
            vars.ykp1_est = obs.ykp1_est;
            vars.p_seq_g_Yk = obs.p_seq_g_Yk;
            vars.gamma_k = obs.gamma_k;
            vars.xkp1_est_f = cell(1, obs.n_filt);
            vars.ykp1_est_f = cell(1, obs.n_filt);
            vars.P_f = cell(1, obs.n_filt);
            for f = 1:obs.n_filt
               vars.xkp1_est_f{f} = obs.filters{f}.xkp1_est;
               vars.ykp1_est_f{f} = obs.filters{f}.ykp1_est;
               vars.P_f{f} = obs.filters{f}.P;
            end
            % Integer variables
            vars.int16.i = obs.i;
            vars.int16.i_next = obs.i_next;

            if strcmp(obs.type, "MKF-SP")
                % Additional variables used by sequence pruning
                % algorithm
                vars.int16.f_main = obs.f_main;
                vars.int16.f_hold = obs.f_hold;
                vars.int16.seq = obs.seq;
            end

        case 'EKF'  % Extended Kalman filters

            % Vars to return
            % TODO: Add dynamic vars

        case 'MEKF'  % Extended Kalman filters

            % Vars to return
            % TODO: Add dynamic vars

        otherwise
            error('Value error: observer type not recognized')
    end

end