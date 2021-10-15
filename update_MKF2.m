function obs = update_MKF2(obs, uk, yk, show_plots)
% obs = update_MKF2(obs) updates the multi-model Kalman
% filter and calculates the estimates of the states and 
% output at the next sample time.
%
% Modification to test lower rate of probability updates.
%
% Arguments:
%   obs : struct containing the multi-model Kalman filter
%       variables (see function mkf_filter).
%   uk : vector (nu, 1) of system inputs at the current 
%       sample time.
%   yk : vector (ny, 1) of system output measurements
%       at the current sample time.
%

    % Debugging option
    if nargin < 4
        show_plots = false;
    end

    % Increment update timer count
    obs.c = obs.c + 1;

    % Do Bayesian update if counter has reached update period
    if obs.c == obs.d

        % Reset counter
        obs.c = 0;

        % Increment sequence index (i should be initialized at 0
        % and 1 <= obs.i <= obs.nf thereafter)
        obs.i = mod(obs.i, obs.nf) + 1;

        % Arrays to store model indicator sequence values
        gamma_k = zeros(obs.n_filt, 1);
        p_gamma_k = nan(obs.n_filt, 1);

        % Bayesian update to conditional probabilities
        for f = 1:obs.n_filt

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
            if ~isscalar(yk_cov)
                yk_cov = triu(yk_cov.',1) + tril(yk_cov);
            end

            % Save for debugging purposes
            obs.filters{f}.yk_cov = yk_cov;

            % Calculate normal probability density (multivariate)
            obs.p_yk_g_seq_Ykm1(f) = mvnpdf(yk, yk_est, yk_cov);

            % Display pdf and yk
            if show_plots && f <= 7
                figure(50+f); clf
                %s3 = 3*sqrt(yk_cov);  % no. of std. dev. to show
                s3 = 1;  % or specify value
                x = linspace(yk_est - s3, yk_est + s3, 101);
                y = normpdf(x, yk_est, sqrt(diag(yk_cov)));
                plot(x, y); hold on
                p_yk_est = normpdf(0, 0, sqrt(diag(yk_cov)));
                stem(yk_est, p_yk_est)
                p_yk = normpdf(yk, yk_est, sqrt(diag(yk_cov)));
                plot(yk, p_yk, 'ok', 'MarkerFaceColor', 'k')
                xlim([yk_est-s3 yk_est+s3])
                ylim([0 5])  % Specify max prob.
                grid on
                title(sprintf('Filter %d',f))
                legend('$p(y(k))$', '$\hat{y}(k)$', '$y_m(k)$','Interpreter','Latex')
                set(gcf,'Position',[f*250-150 50 250 150])
            end

            % Model indicator value from previous sample time
            gamma_km1 = gamma_k(f);

            % Update model indicator value gamma(k) with the
            % current value from the filter's sequence
            gamma_k(f) = obs.seq{f}(:, obs.i);

            % Compute Pr(gamma(k)) based on Markov transition
            % probabilities
            p_gamma_k(f) = prob_gamma(gamma_k(f), obs.T(gamma_km1+1, :)');

            % Select filter system model based on current
            % model indicator value
            ind = gamma_k(f) + 1;
            obs.filters{f}.A = obs.A{ind};
            obs.filters{f}.B = obs.B{ind};
            obs.filters{f}.C = obs.C{ind};
            obs.filters{f}.D = obs.D{ind};
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
    
    end

    % Arrays to collect estimates from each filter
    Xkf_est = zeros(obs.n_filt, obs.n);
    Ykf_est = zeros(obs.n_filt, obs.ny);

    % Update all filter estimates
    for f = 1:obs.n_filt

        % Update observer estimates and covariance matrices
        obs.filters{f} = update_KF(obs.filters{f}, uk, yk);
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