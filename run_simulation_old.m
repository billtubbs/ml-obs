function sim_out = run_simulation_old(Ts, data, u_meas, observers)
% DEPRECATED: Use run_simulation_obs instead.
% sim_out = run_simulation(Ts, data, u_meas, observers)
% Runs a simulation of a set of pre-defined observers fed with
% input data U and output data Y.
%
% Arguments:
%   Ts : sample time.
%   data : table of data including measured input signals named
%       'U', measured output signals named 'Y'. May also contain
%       data for other variables to be used by the observers
%       or to be inlcuded in the simulation outputs.
%   observers : cell array of observer structs.
%
% Returns:
%   sim_out : struct with the following fields: 
%       data : table containing simulation results.
%       if one or more multi-model filters are included:
%       MKF_obs : array of indices indicating which 
%           observers are multi-model observers.
%       MKF_i : filter sequence index.
%       MKF_X_est : array of estimates of each Kalman filter.
%       MKF_p_seq_g_Yk : conditional probabilities.
%

    % Extract data needed for simulation
    nT = size(data, 1) - 1;
    U = data{:, 'U'};
    Y_m = data{:, 'Y_m'};

    % Labels
    data_labels = data.Properties.VariableNames;
    other_data_labels = setdiff(data_labels, {'U', 'Y_m'});

    % Simulation setup
    k = (0:nT)';
    t = Ts*k;

    % Prepare arrays for storing simulation variables
    n_obs = numel(observers);
    X_est = cell(1, n_obs);
    Y_est = cell(1, n_obs);
    for i = 1:n_obs
        n = size(observers{i}.xkp1_est, 1);
        ny = size(observers{i}.ykp1_est, 1);
        X_est{i} = nan(nT+1, n);
        Y_est{i} = nan(nT+1, ny);
    end

    % Find and count number of MKF filters
    n_obs_mkf = 0;
    obs_mkf = nan(1, n_obs_mkf);
    for i = 1:n_obs
        if startsWith(observers{i}.label, "MKF") ...
                || startsWith(observers{i}.label, "AFMM")
            n_obs_mkf = n_obs_mkf + 1;
            obs_mkf(n_obs_mkf) = i;
        end
    end

    % Cell arrays to store data on each MKF observer
    MKF_i = cell(1, n_obs_mkf);
    MKF_X_est = cell(1, n_obs_mkf);
    MKF_Y_est = cell(1, n_obs_mkf);
    MKF_p_seq_g_Yk = cell(1, n_obs_mkf);
    for f = 1:n_obs_mkf
        obs = observers{obs_mkf(f)};
        MKF_i{f} = zeros(size(k));
        MKF_X_est{f} = nan(size(k, 1), obs.n_filt*n);
        MKF_p_seq_g_Yk{f} = nan(size(k, 1), obs.n_filt);
    end

    % Choose debugging point
    k_stop = -1;  % no debugging
    %k_stop = 60;
    show_plots = false;

    % Start simulation at k = 0 (i = 1)
    for ki = 0:nT
        
        i = ki + 1;  % array index starts at 1
        
        % Use this for debugging - set k_stop above
        if ki == k_stop
            show_plots = true;
            k_stop = ki + 1;
        else
            show_plots = false;
        end

        % Get current process inputs
        uk_m = U(i, :)';

        % Get current process outputs
        yk_m = Y_m(i, :)';

        % Record observer estimates from previous timestep
        for j = 1:n_obs
            X_est{j}(i, :) = observers{j}.xkp1_est';
            Y_est{j}(i, :) = observers{j}.ykp1_est';
        end

        % Update observers using measurements uk_m and yk_m
        for j = 1:n_obs
            obs = observers{j};  % makes a copy

            switch obs.type, 

                case {'KF', 'LB'}
                    % Kalman filter or Luenberger observer

                    % Update observer estimates
                    obs = update_KF(obs, uk_m, yk_m);

                    case 'SKF'
                    % Scheduled Kalman filters

                    % Provide actual shock occurence
                    gamma_k = data{i, 'gamma'};

                    % Update observer estimates
                    obs = update_SKF(obs, uk_m, yk_m, gamma_k);

                case {'MKF', 'MKF_RODD'}
                    % Multi-model Kalman filters

                    % Update observer estimates
                    obs = update_MKF(obs, uk_m, yk_m, show_plots);

                    % Save simulation data for plotting later
                    f_mkf = find(obs_mkf == j);
                    MKF_i{f_mkf}(i) = obs.i(1);
                    for f = 1:obs.n_filt
                        MKF_X_est{f_mkf}(i, n*(f-1)+1:n*f) = ...
                            obs.filters{f}.xkp1_est';
                        MKF_Y_est{f_mkf}(i, ny*(f-1)+1:ny*f) = ...
                            obs.filters{f}.ykp1_est';
                    end
                    MKF_p_seq_g_Yk{f_mkf}(i, :) = obs.p_seq_g_Yk';

                case 'MKF_AFMM'
                    % Adaptive forgetting multi-model observers

                    % Update observer estimates
                    obs = update_AFMM(obs, uk_m, yk_m);

                    % Save simulation data for plotting later
                    f_mkf = find(obs_mkf == j);
                    MKF_i{f_mkf}(i) = obs.i(1);
                    for f = 1:obs.n_filt
                        MKF_X_est{f_mkf}(i, n*(f-1)+1:n*f) = ...
                            obs.filters{f}.xkp1_est';
                        MKF_Y_est{f_mkf}(i, ny*(f-1)+1:ny*f) = ...
                            obs.filters{f}.ykp1_est';
                    end
                    MKF_p_seq_g_Yk{f_mkf}(i, :) = obs.p_seq_g_Yk';

                case 'EKF'
                    % Extended Kalman filters

                    % Update observer estimates
                    obs = update_EKF(obs, yk_m, uk_m, Ts);

                case 'MEKF'
                    % Multi-model extended Kalman filters

                    % Update observer estimates
                    obs = update_MEKF(obs, yk_m, uk_m, Ts);

                otherwise
                    error('Value error: observer type not recognized')

            end

            % save changes
            observers{j} = obs;

        end
        
        % Use this for step-by-step debugging
        if show_plots
            fprintf("%d: t = %g\n", ki, t(i));
            disp("k_stop reached");  % debug here
            show_plots = false;
        end
        
    end

    % Record final observer estimates
    for j = 1:n_obs
        X_est{j}(i, :) = observers{j}.xkp1_est';
        Y_est{j}(i, :) = observers{j}.ykp1_est';
    end

    % Calculate output estimation errors
    E_obs = repmat(Y_m, 1, n_obs) - cell2mat(Y_est);
    
    % Convert estimates to tables
    X_est = table(X_est{:}, 'VariableNames', ...
        vector_element_labels('X_est', '', n_obs));
    Y_est = table(Y_est{:}, 'VariableNames', ...
        vector_element_labels('Y_est', '', n_obs));

    % Main simulation results table
    sim_out.data = [table(k,t,U,Y_m) X_est Y_est table(E_obs) ...
        data(:, other_data_labels)];

    % Data on MKF filters
    sim_out.MKF_obs = obs_mkf;
    sim_out.MKF_i = MKF_i;
    sim_out.MKF_X_est = MKF_X_est;
    sim_out.MKF_Y_est = MKF_Y_est;
    sim_out.MKF_p_seq_g_Yk = MKF_p_seq_g_Yk;
    sim_out.observers = observers;

end