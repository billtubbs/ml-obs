function sim_out = run_simulation_obs(Ts, nT, A, B, C, U, alpha, ...
    Pd, V, W, x0, u_meas, y_meas, observers)
% sim_out = run_simulation_obs(Ts, nT, A, B, C, U, alpha, ...
%     V, W, x0, u_meas, y_meas, observers)
% Runs a simulation of the system A, B, C, with inputs U, 
% measurement noise V, and the cell array of observers.
%
% See main script rod_obs_test2.m for details.
%
% Arguments:
%   Ts : sample time.
%   nT : number of sample periods to simulate.
%   n : number of system states.
%   A, B, C : system matrices 
%   U : input signals (incuding unmeasured inputs).
%   alpha : randomly-occuring shock signal.
%   Pd : disturbance signals at input to plant.
%   V : measurement noise signals.
%   W : process noise signals.
%   x0 : Initial values of process states at t=0.
%   u_meas : Boolean vector indicating which inputs are measured
%   y_meas : Boolean vector indicating which outputs are measured
%   observers : cell array of observer structs.
%
% Returns:
%   simout : struct with the following fields: data, and  
%       MKF_i, MKF_X_est, and MKF_p_seq_g_Yk if one or more
%       multi-model filters are included.
%

    assert(size(U, 1) == nT+1);
    assert(size(alpha, 1) == nT+1);
    assert(size(Pd, 1) == nT+1);
    assert(size(V, 1) == nT+1);
    assert(size(W, 1) == nT+1);

    % Dimensions
    n = size(A, 1);
    nu = size(B, 2);
    ny = size(C, 1);
    n_dist_u = sum(~u_meas);

    % Simulation setup
    k_start = -1;  % only needed for high order systems with delays
    k = (k_start:nT)';
    t = Ts*k;

    % Prepare arrays for storing simulation variables
    U = [zeros(-k_start, nu); U];
    alpha = [zeros(-k_start, n_dist_u); alpha];
    V = [zeros(-k_start, ny); V];
    W = [zeros(-k_start, n); W];
    X = zeros(size(k, 1), n);
    Y = zeros(size(k, 1), ny);
    Y_m = zeros(size(Y));
    Pd = [zeros(-k_start, size(Pd, 2)); Pd];
    n_obs = numel(observers);
    X_est = nan(size(k,1), n_obs*n);
    Y_est = nan(size(k,1), n_obs*ny);
    E_obs = zeros(size(k,1), n_obs*ny);

    % Find and count number of MKF filters
    n_obs_mkf = 0;
    observers_mkf = nan(1, n_obs_mkf);
    for i = 1:n_obs
        if startsWith(observers{i}.label, "MKF")
            n_obs_mkf = n_obs_mkf + 1;
            observers_mkf(n_obs_mkf) = i;
        end
    end

    % Cell arrays to store data on each MKF observer
    MKF_i = cell(1, n_obs_mkf);
    MKF_X_est = cell(1, n_obs_mkf);
    MKF_p_seq_g_Yk = cell(1, n_obs_mkf);
    for f = 1:n_obs_mkf
        obs = observers{observers_mkf(f)};
        MKF_i{f} = zeros(size(k));
        MKF_X_est{f} = nan(size(k, 1), obs.n_filt*n);
        MKF_p_seq_g_Yk{f} = nan(size(k, 1), obs.n_filt);
    end

    % Initialize process and observers
    xk = x0;
    for j = 1:n_obs
        observers{j}.xkp1_est = zeros(n,1);
        observers{j}.ykp1_est = zeros(ny,1);
        observers{j}.error = zeros(ny,1);
    end

    % Start simulation at k = 0 (i = 1)
    for ki = 0:nT
        
        i = find(ki == k);
        
        % Calculate current process outputs
        yk = C*xk;
        yk_m = yk + V(i);

        % Record system states and outputs
        X(i, :) = xk';
        Y(i, :) = yk';
        Y_m(i, :) = yk_m';

        % Record observer estimates
        for j = 1:n_obs
            X_est(i, n*(j-1)+1:n*j) = observers{j}.xkp1_est';
            Y_est(i, ny*(j-1)+1:ny*j) = observers{j}.ykp1_est';
        end
        
        % Get current input
        uk = U(i, :)';

        % Input vector for observers with unmeasured inputs
        % set to zero
        uk_m = uk;
        uk_m(~u_meas) = 0;

        % Update observers using measurements uk_m and yk_m
        for j = 1:n_obs

            obs = observers{j};  % Note: makes a copy
            switch obs.label

                case 'none'  % no observer
                    obs.xkp1_est = xk;

                case {'KF', 'KF1', 'KF2', 'KF3'}  % standard Kalman filters

                    % Update observer gains and covariance matrix
                    obs = update_KF(obs, uk_m, yk_m);

                case {'LB', 'LB1', 'LB2', 'KFSS', 'KFSS1', 'KFSS2'}  
                    % steady-state observers

                    % Update observer estimates
                    obs = update_KF(obs, uk_m, yk_m);

                case {'SKF', 'SKF1', 'SKF2'}  % Scheduled Kalman filters

                    % Set process noise covariance matrix Q based on
                    % actual shock occurence
                    a = alpha(i, :);
                    n_dist = size(a, 2);
                    x_var = diag(obs.Q0);
                    x_var(~u_meas) = obs.sigma_wp(sub2ind(size(obs.sigma_wp), ...
                        1:n_dist, a+1)).^2;
                    obs.Q = diag(x_var);
                    
                    % Update observer gains and covariance matrix
                    obs = update_KF(obs, uk_m, yk_m);

                case {'MKF', 'MKF1', 'MKF2', 'MKF3', 'MKF1m', 'MKF2m'}  
                    % multi-model Kalman filters

                    % Update observer estimates
                    obs = update_MKF(obs, uk_m, yk_m);

                    % Save simulation data for plotting later
                    f_mkf = find(observers_mkf == j);
                    MKF_i{f_mkf}(i) = obs.i;
                    for f = 1:obs.n_filt
                        MKF_X_est{f_mkf}(i, n*(f-1)+1:n*f) = ...
                            obs.filters{f}.xkp1_est';
                    end
                    MKF_p_seq_g_Yk{f_mkf}(i, :) = obs.p_seq_g_Yk';

                otherwise
                    error('Value error: observer type not recognized')

            end

            % save changes
            observers{j} = obs;

        end

        % Calculate process states in next timestep
        xk = A*xk + B*uk + W(i, :)';

    end

    % Record final observer estimates
    i = find(k == nT);
    for j = 1:n_obs
        X_est(i, n*(j-1)+1:n*j) = observers{j}.xkp1_est';
        Y_est(i, ny*(j-1)+1:ny*j) = observers{j}.ykp1_est';
    end
    
    % Calculate output estimation errors
    E_obs = repmat(Y, 1, n_obs) - Y_est;
    
    % Main simulation results table
    sim_out.data = table(k,t,U,alpha,V,W,Pd,X,X_est,Y,Y_m,Y_est,E_obs);

    % Remove leading time periods (k < 0)
    sim_out.data = sim_out.data(k >= 0, :);
    for f = 1:n_obs_mkf
        MKF_i{f} = MKF_i{f}(k >= 0);
        MKF_p_seq_g_Yk{f} = MKF_p_seq_g_Yk{f}(k >= 0, :);
        MKF_X_est{f} = MKF_X_est{f}(k >= 0, :);
    end

    % Data on MKF filters
    sim_out.MKF_i = MKF_i;
    sim_out.MKF_X_est = MKF_X_est;
    sim_out.MKF_p_seq_g_Yk = MKF_p_seq_g_Yk;
    sim_out.observers = observers;

end