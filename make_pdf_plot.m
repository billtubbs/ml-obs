function make_pdf_plot(y_mean, y_cov, yk, labels, x_label, sd, y_lim)
% Display normal probability density and a data point.
%
% Arguments:
%   y_mean : double
%     Mean value
%   y_cov : double matrix
%     Covariance
%   yk : double
%     Data point
%   labels : char array
%     Labels to put in legend.
%   x_label : string or char
%     Label for x-axis.
%   sd : (optional, default 3)
%     No. of std. deviations to show.
%   y_lim : (1, 2) row vector (optional)
%     Lower and upper y-axes range.
%

    if nargin < 7
        y_lim = [-inf inf];
    end
    if nargin < 6
        sd = 3;
    end
    if nargin < 5
        x_label = "$y(k),\hat{y}(k),y_M(k)$";
    end
    if nargin < 4
        labels = {'$p(y(k))$', '$\mathrm{Pr}(\hat{y}(k))$', '$\mathrm{Pr}(y_M(k))$'};
    end
    y_min_max = sd * sqrt(y_cov);
    x = linspace(y_mean - y_min_max, y_mean + y_min_max, 101);
    y = normpdf(x, y_mean, sqrt(diag(y_cov)));
    plot(x, y); hold on
    p_yk_est = normpdf(0, 0, sqrt(diag(y_cov)));
    stem(y_mean, p_yk_est)
    p_yk = normpdf(yk, y_mean, sqrt(diag(y_cov)));
    plot(yk, p_yk, 'ok', 'MarkerFaceColor', 'k')
    xlim([y_mean-y_min_max y_mean+y_min_max])
    ylim(y_lim)
    xlabel(x_label,'Interpreter','latex')
    ylabel('Probability','Interpreter','latex')
    grid on
    legend(labels,'Interpreter','latex')

end