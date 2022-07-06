function make_MKF_pdf_plot(obs, f, yk, y_lim)
% Display probability density of the output estimate
% of filter f of the mutli-model observer compared to
% the current data point.
%

    if nargin < 4
        y_lim = [-inf inf];
    end

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

    x_label = "$y(k),\hat{y}(k),y_M(k)$";
    labels = {'$p(y(k))$', '$\hat{y}(k)$', '$y_m(k)$'};
    sd = 3;
    make_pdf_plot(yk_est, yk_cov, yk, labels, x_label, sd, y_lim)

%     %s3 = 3*sqrt(yk_cov);  % no. of std. dev. to show
%     s3 = 1;  % or specify value
%     x = linspace(yk_est - s3, yk_est + s3, 101);
%     y = normpdf(x, yk_est, sqrt(diag(yk_cov)));
%     plot(x, y); hold on
%     p_yk_est = normpdf(0, 0, sqrt(diag(yk_cov)));
%     stem(yk_est, p_yk_est)
%     p_yk = normpdf(yk, yk_est, sqrt(diag(yk_cov)));
%     plot(yk, p_yk, 'ok', 'MarkerFaceColor', 'k')
%     xlim([yk_est-s3 yk_est+s3])
%     ylim([0 5])  % Specify max prob.
%     grid on
%     title(sprintf('Filter %d',f))
%     legend('$p(y(k))$', '$\hat{y}(k)$', '$y_m(k)$','Interpreter','Latex')

end