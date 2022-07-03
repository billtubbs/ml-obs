% Function to run simulations of MKF_SF and MKF_SP
% observers for unit testing.

function [obs, sim_results] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m, ...
    obs)

    k = (0:nT)';
    t = Ts*k;
    X_est = nan(nT+1,n);
    Y_est = nan(nT+1,ny);
    E_obs = nan(nT+1,ny);
    
    % Arrays to store observer variables
    switch obs.type
        case {'KF', 'KFSS'}
            K_obs = cell(nT+1, 1);
        case {'MKF', 'MKF_SF', 'MKF_SF95', 'MKF_SP'}
            n_filt = obs.n_filt;
            MKF_p_seq_g_Yk = nan(nT+1, n_filt);
            K_obs_j = cell(nT+1, n_filt);
            trP_obs_j = cell(nT+1, n_filt);
    end
    trP_obs = cell(nT+1, 1);

    % Start simulation at k = 0
    for i = 1:nT+1

        % For debugging:
        %fprintf("t = %f\n", t(i));

        % Process measurements
        uk_m = U_m(i,:)';
        yk_m = Y_m(i,:)';

        % Record observer estimates and output errors
        X_est(i, :) = obs.xkp1_est';
        Y_est(i, :) = obs.ykp1_est';
        trP_obs{i, 1} = trace(obs.P);
        E_obs(i, :) = yk_m' - obs.ykp1_est';

        switch obs.type

            case {'KF', 'SKF'}

                % Record filter gain and covariance matrix
                K_obs{i, 1} = obs.K';

            case {'MKF', 'MKF_SF', 'MKF_SF95'}

                % Record filter gains and covariance matrices of
                % each model filter
                for j = 1:obs.n_filt
                    K_obs_j{i, j} = obs.filters{j}.K';
                    trP_obs_j{i, j} = trace(obs.filters{j}.P);
                end

                % Record filter conditional probabilities
                MKF_p_seq_g_Yk(i, :) = obs.p_seq_g_Yk';

            case {'MKF_SP'}

                % Record filter gains and covariance matrices of
                % each model filter
                for j = 1:obs.n_filt
                    K_obs_j{i, j} = obs.filters{j}.K';
                    trP_obs_j{i, j} = trace(obs.filters{j}.P);
                end

                % Record filter conditional probabilities
                MKF_p_seq_g_Yk(i, :) = obs.p_seq_g_Yk';

                % Record filter arrangement
                MKF_SP_f_main(i, :) = obs.f_main;
                MKF_SP_f_hold(i, :) = obs.f_hold;

            otherwise
                error('Observer type not valid')


        end

        % Update observer and estimates
        obs.update(yk_m, uk_m);

    end

    sim_results.t = t;
    sim_results.k = k;
    sim_results.X_est = X_est;
    sim_results.Y_est = Y_est;
    sim_results.E_obs = E_obs;
    sim_results.trP_obs = trP_obs;
    switch obs.type
        case {'KF', 'KFSS'}
            sim_results.K_obs = K_obs;
        case {'MKF', 'MKF_SF', 'MKF_SF95'}
            sim_results.MKF_p_seq_g_Yk = MKF_p_seq_g_Yk;
            sim_results.K_obs_j = K_obs_j;
            sim_results.trP_obs_j = trP_obs_j;
        case 'MKF_SP'
            sim_results.MKF_SP_f_main = MKF_SP_f_main;
            sim_results.MKF_SP_f_hold = MKF_SP_f_hold;
    end

end