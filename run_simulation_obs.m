function [Xk_est,Yk_est,DiagPk,MKF_vars] = ...
    run_simulation_obs(Ym,U,Gamma,seq,observers,f_mkf)
% [Xk_est,Yk_est,DiagPk,MKF_vars] = ...
%     run_simulation_obs(Ym,U,Gamma,seq,observers,f_mkf)
% Simulate a set of observers on a switching system. Used by the 
% following test scripts:
%
%  - test_MKFObservers_JS.m
%  - test_MKFObservers_2x2.m
%

    nT = size(Ym, 1) - 1;
    ny = size(Ym, 2);
    n_obs = numel(observers);
    n = size(observers{1}.xk_est, 1);

    Xk_est = zeros(nT+1, n*n_obs);
    Yk_est = zeros(nT+1, ny*n_obs);
    DiagPk = zeros(nT+1, n*n_obs);

    if ~isempty(f_mkf)
        obs_mkf = observers{f_mkf};
        MKF_label = obs_mkf.label;
        nh = obs_mkf.nh;
        MKF_Xk_est = cell(nT+1, nh);
        MKF_Yk_est = cell(nT+1, nh);
        MKF_K_obs = cell(nT+1, nh);
        MKF_trP_obs = nan(nT+1, nh);
        MKF_i = nan(nT+1, 1);
        MKF_p_seq_g_Yk = nan(nT+1, nh);
    else
        MKF_Xk_est = {};
        MKF_Yk_est = {};
        MKF_K_obs = {};
        MKF_trP_obs = nan;
        MKF_i = nan;
        MKF_p_seq_g_Yk = nan;
    end

    if iscell(seq)
        seq = cell2mat(seq);  % makes it easier to index
    end

    for i = 1:nT+1

        yk = Ym(i, :)';
        uk = U(i, :)';

        % Update observers
        for f = 1:n_obs
            obs = observers{f};
            switch obs.type
                case "SKF"
                    rk = Gamma(i) + 1;
                    obs.update(yk, uk, rk);
                case {"MKF", "MKF_F"}
                    rk = seq(:, i);
                    obs.update(yk, uk, rk);
                otherwise
                    obs.update(yk, uk);
            end
            if f == f_mkf
                if isprop(obs, "i")
                    MKF_i(i, :) = obs.i;
                end
                MKF_p_seq_g_Yk(i, :) = obs.p_seq_g_Yk';
                for j = 1:obs.nh
                    MKF_Xk_est{i, j} = obs.filters.Xk_est(:, :, j);
                    MKF_Yk_est{i, j} = obs.filters.Yk_est(:, :, j);
                    MKF_K_obs{i, j} = obs.filters.Kf(:, :, j);
                    MKF_trP_obs(i, j) = trace(obs.filters.Pkp1(:, :, j));
                end
            end
            xk_est(1, (f-1)*n+1:f*n) = obs.xk_est';
            yk_est(1, (f-1)*ny+1:f*ny) = obs.yk_est';
            diagP(1, (f-1)*n+1:f*n) = diag(obs.Pk)';
        end

        % Record observer estimates
        Xk_est(i, :) = xk_est;
        Yk_est(i, :) = yk_est;
        DiagPk(i, :) = diagP;

        % Store MKF results in struct
        MKF_vars.label = MKF_label;
        MKF_vars.Xk_est = MKF_Xk_est;
        MKF_vars.Yk_est = MKF_Yk_est;
        MKF_vars.K_obs = MKF_K_obs;
        MKF_vars.trP_obs= MKF_trP_obs;
        MKF_vars.i = MKF_i;
        MKF_vars.p_seq_g_Yk= MKF_p_seq_g_Yk;

    end
end