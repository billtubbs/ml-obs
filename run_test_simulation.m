function [obs, sim_results] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m, ...
    obs,show_plots)
% Function to run test simulations of MKF_SF and MKF_SP
% observers. Records values of various internal variables
% for unit testing, analysis, and debugging.
%

    if nargin < 8
        show_plots = false;
    end
    k = (0:nT)';
    t = Ts*k;
    X_est = nan(nT+1,n);
    Y_est = nan(nT+1,ny);
    E_obs = nan(nT+1,ny);
    
    % Arrays to store observer variables
    switch obs.type
        case {'KF', 'KFSS', 'SKF'}
            K_obs = cell(nT+1, 1);
        case {'MKF', 'MKF_SF', 'MKF_SF95', 'MKF_SP'}
            n_filt = obs.n_filt;
            MKF_p_seq_g_Yk = nan(nT+1, n_filt);
            K_obs_f = cell(nT+1, n_filt);
            trP_obs_f = cell(nT+1, n_filt);
            MKF_X_est_f = cell(nT+1, n_filt);
        case {'MKF_SP'}
            n_filt = obs.n_filt;
            MKF_p_seq_g_Yk = nan(nT+1, n_filt);
            K_obs_f = cell(nT+1, n_filt);
            trP_obs_f = cell(nT+1, n_filt);
            MKF_X_est_f = cell(nT+1, n_filt);
            MKF_SP_f_main(i, :) = int16(zeros(nT+1, obs.f_main));
            MKF_SP_f_hold(i, :) = int16(zeros(nT+1, obs.f_hold));
        otherwise
            error("Observer type not recognized")
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

                % Record filter gains, trace of covariance matrices
                % and state estimates of each model filter
                for j = 1:obs.n_filt
                    K_obs_f{i, j} = reshape(obs.filters{j}.K, 1, []);
                    trP_obs_f{i, j} = trace(obs.filters{j}.P);
                    MKF_X_est_f{i, j} = obs.filters{j}.xkp1_est';
                end

                % Record filter conditional probabilities
                MKF_p_seq_g_Yk(i, :) = obs.p_seq_g_Yk';

            case {'MKF_SP'}

                % Record filter gains, trace of covariance matrices
                % and state estimates of each model filter
                for j = 1:obs.n_filt
                    K_obs_f{i, j} = obs.filters{j}.K';
                    trP_obs_f{i, j} = trace(obs.filters{j}.P);
                    MKF_X_est_f{i, j} = obs.filters{j}.xkp1_est';
                end

                % Record filter conditional probabilities
                MKF_p_seq_g_Yk(i, :) = obs.p_seq_g_Yk';

                % Record filter arrangement
                MKF_SP_f_main(i, :) = obs.f_main;
                MKF_SP_f_hold(i, :) = obs.f_hold;

            otherwise
                error('Observer type not valid')

        end

        if show_plots && t(i) >= 34.5
            for f = 1:min(obs.n_filt, 8)
                figure(50+f); clf
                make_MKF_pdf_plot(obs, f, yk_m, [-2 2], 5)
                title(sprintf('Filter %d',f))
                legend('$p(y(k))$', '$\hat{y}(k)$', '$y_m(k)$','Interpreter','Latex')
                set(gcf, 'Position', [f*250-150 100 250 150])
            end
        end

        % Update observer and estimates
        obs.update(yk_m, uk_m);

        if show_plots && t(i) >= 34.5
            % Plot gamma_k, p_yk_g_seq_Ykm1 and 
            figure(50)
            subplot(4,1,1)
            bar(1:obs.n_filt, obs.p_gamma_k)
            ylim([0, 1])
            title('Shock probabilities $Pr(\gamma(k) \mid \gamma(k-1))$','Interpreter','Latex')
            subplot(4,1,2)
            bar(1:obs.n_filt, obs.p_seq_g_Ykm1)
            title('$p(\Gamma(k) \mid Y(k-1))$','Interpreter','Latex')
            ylim([0, 1])
            subplot(4,1,3)
            bar(1:obs.n_filt, obs.p_yk_g_seq_Ykm1)
            title('$p(y_M(k) \mid \Gamma(k),Y(k-1))$','Interpreter','Latex')
            subplot(4,1,4)
            bar(1:obs.n_filt, obs.p_seq_g_Yk)
            xlabel('Filter')
            title('$p(\Gamma(k) \mid Y(k))$','Interpreter','Latex')
            fprintf("t(%d): %g\n", k(i), t(i))
            if t(i) >= 34.5
                txt = input("Paused. Enter 'd' to debug, 'q' to quit: ", "s");
                if strcmpi(txt,"D")
                    dbstop if true
                elseif strcmpi(txt,"Q")
                    return
                end
            end
        end

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
            sim_results.K_obs_f = K_obs_f;
            sim_results.trP_obs_f = trP_obs_f;
            sim_results.MKF_X_est_f = MKF_X_est_f;
        case 'MKF_SP'
            sim_results.MKF_p_seq_g_Yk = MKF_p_seq_g_Yk;
            sim_results.K_obs_f = K_obs_f;
            sim_results.trP_obs_f = trP_obs_f;
            sim_results.MKF_X_est_f = MKF_X_est_f;
            sim_results.MKF_SP_f_main = MKF_SP_f_main;
            sim_results.MKF_SP_f_hold = MKF_SP_f_hold;
    end

end