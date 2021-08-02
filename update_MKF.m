function obs = update_MKF(obs, uk, yk)
% obs = update_MKF(obs) updates the multi-model Kalman
% filter and calcultes the estimates of the states and 
% output in current time period.
%
% Arguments:
%   obs : struct containing the multi-model Kalman filter
%       variables (see function mkf_filter).
%   uk : vector (nu, 1) of system inputs in current 
%       sample period.
%   yk : vector (ny, 1) of system output measurements
%       in current sample period.
%

    % Debugging option
    show_plots = false;
    
    % Update sequence index (i should be initialized at 0
    % and 1 <= obs.i <= obs.nf)
    obs.i = mod(obs.i, obs.nf) + 1;

    % Arrays to collect estimates from each filter
    Xkf_est = zeros(obs.n_filt, obs.n);
    Ykf_est = zeros(obs.n_filt, obs.ny);
    gamma_k = zeros(obs.n_filt, 1);
    p_gamma_k = zeros(obs.n_filt, 1);

    % Update all filters
    for f = 1:obs.n_filt
        
        % TODO: Is the timing/order here correct?

        % Compute posterior probability density of y(k)
        % using posterior PDF (normal distribution) and
        % estimates computed in previous timestep

        % Get y_est(k/k-1) estimated in previous time step
        yk_est = obs.filters{f}.ykp1_est;

        % Calculate covariance of the output estimation errors
        P = obs.filters{f}.P;
        C = obs.filters{f}.C;
        yk_cov = C*P*C' + obs.filters{f}.R;

        % Make sure covariance matrix is symmetric
        yk_cov = triu(yk_cov.',1) + tril(yk_cov);

        % Save for debugging purposes
        obs.filters{f}.yk_cov = yk_cov;

        % Calculate normal probability density (multivariate)
        obs.p_yk_g_seq_Ykm1(f) = mvnpdf(yk, yk_est, yk_cov);

        % Display pdf and yk
        if show_plots && f <= 7
            figure(50+f); clf
            %s3 = 3*sqrt(yk_cov);
            s3 = 0.5;
            x = linspace(yk_est - s3, yk_est + s3, 101);
            y = normpdf(x, yk_est, sqrt(diag(yk_cov)));
            plot(x, y); hold on
            p_yk_est = normpdf(0, 0, sqrt(diag(yk_cov)));
            stem(yk_est, p_yk_est)
            p_yk = normpdf(yk, yk_est, sqrt(diag(yk_cov)));
            plot(yk, p_yk, 'ok', 'MarkerFaceColor', 'k')
            xlim([yk_est-s3 yk_est+s3])
            ylim([0 10])
            grid on
            title(sprintf('Filter %d',f))
            legend('$p(y(k))$', '$\hat{y}(k)$', '$y_m(k)$','Interpreter','Latex')
            set(gcf,'Position',[f*250-150 50 250 250])
        end

        % Current model indicator values from each
        % filter's sequence
        gamma_k(f) = obs.seq{f}(:, obs.i);

        % Compute Pr(gamma(k)) which is assumed to be an
        % independent random variable
        p_gamma_k(f) = prob_gamma(gamma_k(f), obs.p_gamma);

        % Select filter system model based on current
        % model indicator value
        ind = gamma_k(f) + 1;
        obs.filters{f}.A = obs.A{ind};
        obs.filters{f}.B = obs.B{ind};
        obs.filters{f}.C = obs.C{ind};
        obs.filters{f}.D = obs.D{ind};
        obs.filters{f}.Q = obs.Q{ind};
        obs.filters{f}.R = obs.R{ind};

        % Update observer estimates, gain and covariance matrix
        obs.filters{f} = update_KF(obs.filters{f}, uk, yk);
        assert(~any(isnan(obs.filters{f}.xkp1_est)))

        % Save state and output estimates for next timestep
        Xkf_est(f, :) = obs.filters{f}.xkp1_est';
        Ykf_est(f, :) = obs.filters{f}.ykp1_est';

    end

    assert(~any(isnan(obs.p_yk_g_seq_Ykm1)))
    assert(~all(obs.p_yk_g_seq_Ykm1 == 0))

    % Compute Pr(Gamma(k)|Y(k-1)) in current timestep from
    % previous estimate (Pr(Gamma(k-1)|Y(k-1)))
    p_seq_g_Ykm1 = p_gamma_k .* obs.p_seq_g_Yk;

    % Bayesian update of Pr(Gamma(k)|Y(k))
    cond_pds = obs.p_yk_g_seq_Ykm1 .* p_seq_g_Ykm1;
    obs.p_seq_g_Yk = cond_pds ./ sum(cond_pds);
    
    % Save variables for debugging purposes
    obs.p_gamma_k = p_gamma_k;
    obs.p_seq_g_Ykm1 = p_seq_g_Ykm1;
    
    if show_plots
        % Plot gamma_k, p_yk_g_seq_Ykm1 and 
        figure(50)
        subplot(3,1,1)
        bar(1:obs.n_filt, obs.p_yk_g_seq_Ykm1)
        title('$p(y(k)|\Gamma(k),Y(k-1))$','Interpreter','Latex')
        subplot(3,1,2)
        bar(1:obs.n_filt, p_gamma_k)
        ylim([0, 1])
        title('Shock probabilities $\gamma(k)$','Interpreter','Latex')
        subplot(3,1,3)
        bar(1:obs.n_filt, obs.p_seq_g_Yk)
        xlabel('Filter')
        title('$p(\Gamma(k)|Y(k))$','Interpreter','Latex')
    end

        % Compute multi-model observer state and output estimates
        obs.xkp1_est = sum(Xkf_est .* obs.p_seq_g_Yk, 1)';
        obs.ykp1_est = sum(Ykf_est .* obs.p_seq_g_Yk, 1)';
        assert(~any(isnan(obs.xkp1_est)))
        assert(~any(isnan(obs.ykp1_est)))

end