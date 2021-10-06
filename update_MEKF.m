function obs = update_MEKF(obs, yk, uk)
% obs = update_MEKF(obs) updates the multi-model extended 
% Kalman filter and calculates the estimates of the states
% and output at the next sample time.
%
% Arguments:
%   obs : struct containing the multi-model Kalman filter
%       variables (see function mkf_filter).
%   uk : vector (nu, 1) of system inputs at the current 
%       sample time.
%   yk : vector (ny, 1) of system output measurements
%       at the current sample time.
%

    % increment sequence index (i should be initialized at 0
    % and 1 <= obs.i <= obs.nf thereafter)
    obs.i = mod(obs.i, obs.nf) + 1;

    % Arrays to collect estimates from each filter
    Xkf_est = zeros(obs.n_filt, obs.n);
    Ykf_est = zeros(obs.n_filt, obs.ny);
    gamma_k = zeros(obs.n_filt, 1);
    p_gamma_k = zeros(obs.n_filt, 1);

    % Update all filters
    for f = 1:obs.n_filt

        % Compute posterior probability density of y(k)
        % using posterior PDF (normal distribution) and
        % estimates computed in previous timestep

        % Get y_est(k/k-1) estimated in previous time step
        yk_est = obs.filters{f}.ykp1_est;

        % Model indicator value from previous sample time
        gamma_km1 = gamma_k(f);

        % Update model indicator value gamma(k) with the
        % current value from the filter's sequence
        gamma_k(f) = obs.seq{f}(:, obs.i);

        % Model index
        ind = gamma_k(f) + 1;

        % Calculate covariance of the output estimation errors
        P = obs.filters{f}.P;
        % Calculate Jacobian of measurement function linearized at
        % current state estimates.
        params = obs.params{ind};  % current model parameters
        H = obs.filters{f}.dhdx(obs.filters{f}.xkp1_est, uk, params);
        yk_cov = H*P*H' + obs.filters{f}.R;

        % Make sure covariance matrix is symmetric
        yk_cov = triu(yk_cov.',1) + tril(yk_cov);

        % Save for debugging purposes
        obs.filters{f}.yk_cov = yk_cov;

        % Calculate normal probability density (multivariate)
        obs.p_yk_g_seq_Ykm1(f) = mvnpdf(yk, yk_est, yk_cov);

        % Compute Pr(gamma(k)) based on Markov transition
        % probabilities
        p_gamma_k(f) = prob_gamma(gamma_k(f), obs.T(gamma_km1+1, :)');

        % Select filter system model based on index value
        obs.filters{f}.f = obs.f{ind};
        obs.filters{f}.h = obs.h{ind};
        obs.filters{f}.dfdx = obs.dfdx{ind};
        obs.filters{f}.dhdx = obs.dhdx{ind};
        obs.filters{f}.Q = obs.Q{ind};
        obs.filters{f}.R = obs.R{ind};

        % Update observer estimates, gain and covariance matrix
        obs.filters{f} = update_EKF(obs.filters{f}, yk, uk, params);
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

end