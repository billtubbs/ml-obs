function obs = set_obs_vars(obs, vars)
% obs = set_obs_vars(obs, vars) sets all the variables (i.e.
% time-varying properties) of the observer from the values
% in vars (a struct).
%

% TODO: Write a test script for set_obs_vars and get_obs_vars
    switch(obs.type)
        
        case 'none'  % no observer

            % No variables to set

        case {'KFSS', 'LB'}  % Steady-state filters

            % Set variables
            obs.xkp1_est = vars.xkp1_est;
            obs.ykp1_est = vars.ykp1_est;

        case 'KF'  % Standard Kalman filters

            % Set variables
            obs.xkp1_est = vars.xkp1_est;
            obs.ykp1_est = vars.ykp1_est;
            obs.P = vars.P;

        case 'SKF'  % Scheduled Kalman filters

            % Set variables
            obs.xkp1_est = vars.xkp1_est;
            obs.ykp1_est = vars.ykp1_est;
            obs.P = vars.P;

        case {'MKF','MKF_SF95', 'MKF_SF'}  % Multi-model Kalman filters

            % Set double variables
            obs.xkp1_est = vars.xkp1_est;
            obs.ykp1_est = vars.ykp1_est;
            obs.p_seq_g_Yk = vars.p_seq_g_Yk;
            obs.gamma_k = vars.gamma_k;
            for f = 1:obs.n_filt
               obs.filters{f}.xkp1_est = vars.xkp1_est_f{f};
               obs.filters{f}.ykp1_est = vars.ykp1_est_f{f};
               obs.filters{f}.P = vars.P_f{f};
            end
            % Set integer variables
            obs.i = vars.int16.i;
            obs.i_next = vars.int16.i_next;

        case 'MKF_SP'  % MKF observer with sequence pruning

            % Set double variables
            obs.xkp1_est = vars.xkp1_est;
            obs.ykp1_est = vars.ykp1_est;
            obs.p_seq_g_Yk = vars.p_seq_g_Yk;
            obs.gamma_k = vars.gamma_k;
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
            obs.seq = vars.int16.seq;

        case 'EKF'  % Extended Kalman filters

            % Set variables
            % TODO: Add dynamic vars

        case 'MEKF'  % Extended Kalman filters

            % Set variables
            % TODO: Add dynamic vars

        otherwise
            error('Value error: observer type not recognized')
    end

end