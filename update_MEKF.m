function obs = update_MEKF(obs, yk, varargin)
% obs = update_MEKF(obs, yk, varargin)
% updates the multi-model extended Kalman filter and 
% calculates the estimates of the states and output at
% the next sample time.
%
% Arguments:
%   obs : struct containing the multi-model Kalman filter
%       variables (see function mkf_filter).
%   uk : vector (nu, 1) of system inputs at the current 
%       sample time.
%   yk : vector (ny, 1) of system output measurements
%       at the current sample time.
%

    % Increment sequence index and update counter
    % obs.i(1) is the sequence index (1 <= i(1) <= obs.f)
    % obs.i(2) is the counter for prob. updates (1 <= i(2) <= obs.d)
    % Whenever obs.i(2) exceeds obs.d (the spacing parameter), it is
    % reset to 1, and the sequence index obs.i(1) is incremented.
    obs.i = obs.i_next;
    obs.i_next = [mod(obs.i(1) - 1 + ...
                  idivide(obs.i(2), obs.d), obs.f) + 1, ...
                  mod(obs.i(2), obs.d) + 1];

    % Do Bayesian update if counter will reset to 1 at next
    % sample time.
    if obs.i_next(2) == 1

        % Arrays to store model indicator sequence values
        gamma_k = zeros(obs.n_filt, 1);
        p_gamma_k = nan(obs.n_filt, 1);

        % Bayesian update to conditional probabilities
        for f = 1:obs.n_filt

            %TODO: I think this should be updated each timestep
            %      only model indicator and covariances are changed
            %      every detection interval.
            % Compute posterior probability density of y(k)
            % using posterior PDF (normal distribution) and
            % estimates computed in previous timestep

            % Get y_est(k/k-1) estimated in previous time step
            yk_est = obs.filters{f}.ykp1_est;

            % Model indicator value from previous sample time
            gamma_km1 = gamma_k(f);

            % Update model indicator value gamma(k) with the
            % current value from the filter's sequence
            gamma_k(f) = obs.seq{f}(:, obs.i(1));

            % Compute Pr(gamma(k)) based on Markov transition
            % probability matrix
            p_gamma_k(f) = prob_gamma(gamma_k(f), obs.T(gamma_km1+1, :)');

            % Model index
            ind = gamma_k(f) + 1;

            % Calculate covariance of the output estimation errors
            P = obs.filters{f}.Pkp1;
            % Calculate Jacobian of measurement function linearized at
            % current state estimates.
            params = obs.params{ind};  % current model parameters
            varargin2 = [varargin obs.params{ind}{:}];
            H = obs.filters{f}.dhdx(obs.filters{f}.xkp1_est, varargin2{:});
            yk_cov = H*P*H' + obs.filters{f}.R;

            % Make sure covariance matrix is symmetric
            if ~isscalar(yk_cov)
                yk_cov = triu(yk_cov.',1) + tril(yk_cov);
            end

            % Save for debugging purposes
            obs.filters{f}.yk_cov = yk_cov;

            % Calculate normal probability density (multivariate)
            obs.p_yk_g_seq_Ykm1(f) = mvnpdf(yk, yk_est, yk_cov);

            % Select filter system model based on index value
            obs.filters{f}.state_fcn = obs.state_fcn{ind};
            obs.filters{f}.meas_fcn = obs.meas_fcn{ind};
            obs.filters{f}.dfdx = obs.dfdx{ind};
            obs.filters{f}.dhdx = obs.dhdx{ind};
            obs.filters{f}.params = obs.params{ind};
            obs.filters{f}.Q = obs.Q{ind};
            obs.filters{f}.R = obs.R{ind};

        end

        assert(~any(isnan(obs.p_yk_g_seq_Ykm1)))
        assert(~all(obs.p_yk_g_seq_Ykm1 == 0))

        % Compute Pr(Gamma(k)|Y(k-1)) in current timestep from
        % previous estimate (Pr(Gamma(k-1)|Y(k-1))) and transition
        % probabilities
        p_seq_g_Ykm1 = p_gamma_k .* obs.p_seq_g_Yk;

        % Bayesian update of Pr(Gamma(k)|Y(k))
        cond_pds = obs.p_yk_g_seq_Ykm1 .* p_seq_g_Ykm1;
        obs.p_seq_g_Yk = cond_pds ./ sum(cond_pds);
        % Note: above calculation normalizes p_seq_g_Yk so that
        % assert(abs(sum(obs.p_seq_g_Yk) - 1) < 1e-15) % is always true

        % Save variables for debugging purposes
        obs.p_gamma_k = p_gamma_k;
        obs.p_seq_g_Ykm1 = p_seq_g_Ykm1;

    end

    % Arrays to collect estimates from each filter
    Xkf_est = zeros(obs.n_filt, obs.n);
    Ykf_est = zeros(obs.n_filt, obs.ny);    

    % Update all filters
    for f = 1:obs.n_filt

        % Update observer estimates, gain and covariance matrix
        obs.filters{f} = update_EKF(obs.filters{f}, yk, varargin{:});
        assert(~any(isnan(obs.filters{f}.xkp1_est)))

        % Save state and output estimates for next timestep
        Xkf_est(f, :) = obs.filters{f}.xkp1_est';
        Ykf_est(f, :) = obs.filters{f}.ykp1_est';

    end

    % Compute multi-model observer state and output estimates
    obs.xkp1_est = sum(Xkf_est .* obs.p_seq_g_Yk, 1)';
    obs.ykp1_est = sum(Ykf_est .* obs.p_seq_g_Yk, 1)';
    assert(~any(isnan(obs.xkp1_est)))
    assert(~any(isnan(obs.ykp1_est)))
    % TODO: Calculate multi-model state covariance estimate

end