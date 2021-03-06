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

    % Arrays to store model indicator sequence values
    gamma_k = zeros(obs.n_filt, 1);
    p_gamma_k = nan(obs.n_filt, 1);

    % Arrays to collect estimates from each filter
    Xkf_est = zeros(obs.n_filt, obs.n);
    Ykf_est = zeros(obs.n_filt, obs.ny);    

    % Bayesian update to conditional probabilities
    for f = 1:obs.n_filt

        % Compute posterior probability density of y(k)
        % using posterior PDF (normal distribution) and
        % estimates computed in previous timestep

        % Get y_est(k/k-1) estimated in previous time step
        yk_est = obs.filters{f}.ykp1_est;

        % Index of model used in previous sample time
        ind_km1 = gamma_k(f) + 1;  % MATLAB indexing

        % Calculate covariance of the output estimation errors
        P = obs.filters{f}.Pkp1;
        % Calculate Jacobian of measurement function linearized at
        % current state estimates.
        % TODO: Should the params from past or current timestep be used?
        varargin2 = [varargin obs.params{ind_km1}{:}];
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

        % Update model indicator value gamma(k) with the
        % current value from the filter's sequence
        gamma_k(f) = obs.seq{f}(:, obs.i(1));

        % Model index at current sample time
        ind = gamma_k(f) + 1;  % MATLAB indexing

        % Compute Pr(gamma(k)) based on Markov transition
        % probability matrix
        p_gamma_k(f) = prob_gamma(gamma_k(f), obs.T(ind_km1, :)');

        % Update filter covariances if at start of a detection
        % interval, TODO: Is this the right time/place to do
        % this update?
        if obs.i(2) == 1  % or should this be obs.i_next(2) == 1

            % Select filter system model based on current
            % model indicator value
            obs.filters{f}.state_fcn = obs.state_fcn{ind};
            obs.filters{f}.meas_fcn = obs.meas_fcn{ind};
            obs.filters{f}.dfdx = obs.dfdx{ind};
            obs.filters{f}.dhdx = obs.dhdx{ind};
            obs.filters{f}.params = obs.params{ind};
            obs.filters{f}.Q = obs.Q{ind};
            obs.filters{f}.R = obs.R{ind};

        end

        % Update observer estimates, gain and covariance matrix
        obs.filters{f} = update_EKF(obs.filters{f}, yk, varargin{:});
        assert(~any(isnan(obs.filters{f}.xkp1_est)))

        % Save state and output estimates for next timestep
        Xkf_est(f, :) = obs.filters{f}.xkp1_est';
        Ykf_est(f, :) = obs.filters{f}.ykp1_est';

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

    % Compute multi-model observer state and output estimates
    obs.xkp1_est = sum(Xkf_est .* obs.p_seq_g_Yk, 1)';
    obs.ykp1_est = sum(Ykf_est .* obs.p_seq_g_Yk, 1)';
    assert(~any(isnan(obs.xkp1_est)))
    assert(~any(isnan(obs.ykp1_est)))
    % TODO: Calculate multi-model state covariance estimate

end