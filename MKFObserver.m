% Multi-model Kalman Filter class definition
%
% obs = MKFObserver(models,P0,T,r0,label,x0,p_seq_g_Yk_init,reset)
% Class for simulating a multi-model Kalman filter for state
% estimation of a Markov jump linear system. 
% 
% This is the filtering form of the observer, which 
% produces posterior estimates of the states and outputs 
% at the current time instant given the data at the 
% current time:
%
%  x_hat(k|k) : estimate of states at time k
%  y_hat(k|k) : estimate of outputs at time k
%
% Prior estimates of the states and outputs at the next
% time instant given the data at the current time are
% also calculated:
%
%  x_hat(k+1|k) : estimate of states at time k + 1
%  y_hat(k+1|k) : estimate of outputs at time k + 1
%
% The system model is defined as:
%
%   x(k+1) = A(k) x(k) + B(k) u(k) + w(k)
%     y(k) = C(k) x(k) + v(k)
%
% Note: there is no direct transmission (D = 0).
%
% Arguments:
%   models : (1, nj) cell array of structs
%       Each struct contains the parameters of a linear
%       model of the system dynamics. These include: A, B, 
%       and C for the system matrices, Q and R for the
%       state error covariance and output measurement 
%       noise covariance, and Ts for the sample period.
%   P0 : (n, n) double
%       Initial covariance matrix of the state estimates
%       (same for each filter).
%   T : Transition probabity matrix of the Markov switching
%       process.
%   label : string
%       Arbitrary name to identify the observer.
%   x0 : (n, 1) double (optional, default zeros)
%       Initial state estimates.
%   r0 : (nh, 1) integer (optional, default ones)
%       Integer in the range {1, ..., nj} which indicates
%       the prior system mode at time k = -1.
%   p_seq_g_Yk_init : (optional, default uniform)
%       Initial prior hypothesis probabilities at time k-1.
%       If not specified, default is equal, i.e. uniform,
%       probability assigned to each hypothesis.
%   reset : logical (default, true)
%       If true, the objects reset method is called after
%       initialization (this is mainly intended for use by
%       other objects instantiating an instance without
%       reseting).
%

classdef MKFObserver < matlab.mixin.Copyable
    properties (SetAccess = immutable)
        Ts (1, 1) double {mustBeNonnegative}
        n (1, 1) double {mustBeInteger, mustBeNonnegative}
        nu (1, 1) double {mustBeInteger, mustBeNonnegative}
        ny (1, 1) double {mustBeInteger, mustBeNonnegative}
        nj (1, 1) double {mustBeInteger, mustBeNonnegative}
        nh (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        models (1, :) cell
        P0 double
        T double
        r0 (:, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(r0, 1)}
        label (1, 1) string
        x0 (:, 1) double
        p_seq_g_Yk_init double
        p_seq_g_Ykm1 double
        p_seq_g_Yk double
        p_yk_g_seq_Ykm1 double
        p_rk_g_Ykm1 double
        p_rk_g_rkm1 double
        filters  struct
        xk_est (:, 1) double
        Pk double
        yk_est (:, 1) double
        xkp1_est (:, 1) double
        Pkp1 double
        rk (:, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(rk, 1)}
        rkm1 (:, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(rkm1, 1)}
        type (1, 1) string
    end
    methods
        function obj = MKFObserver(models,P0,T,r0,label,x0, ...
                p_seq_g_Yk_init,reset)
            arguments
                models (1, :) cell
                P0 double
                T double
                r0 (:, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(r0, 1)}
                label (1, 1) string = ""
                x0 = []
                p_seq_g_Yk_init = []
                reset logical = true
            end

            % Get number of system models and check their dimensions
            [nj, n, nu, ny, Ts] = check_models(models);

            % Number of hypotheses to be modelled
            nh = size(r0, 1);

            % Check dimensions of other parameters
            assert(isequal(size(P0), [n n]), "ValueError: size(P0)")
            for j = 1:nj
                assert(isequal(size(models{j}.Q), [n n]))
                assert(isequal(size(models{j}.R), [ny ny]))
            end

            if isempty(p_seq_g_Yk_init)
                % Initial values of prior conditional probabilities at 
                % k = -1. In absence of prior knowledge, assume all 
                % equally likely
                p_seq_g_Yk_init = ones(nh, 1) ./ nh;
            end

            % Determine initial state values
            if isempty(x0)
                x0 = zeros(n, 1);  % default initial states
            else
                assert(isequal(size(x0), [n 1]), "ValueError: x0")
            end

            % Check transition probability matrix
            assert(isequal(size(T), [nj nj]), "ValueError: size(T)")
            assert(all(abs(sum(obj.T, 2) - 1) < 1e-15), "ValueError: T")

            % Create struct to store Kalman filter variables
            obj.filters = struct();
            obj.filters.Xkp1_est = nan(n, 1, obj.nh);
            obj.filters.Pkp1 = nan(n, n, obj.nh);
            obj.filters.Xk_est = nan(obj.n, 1, obj.nh);
            obj.filters.Pk = nan(obj.n, obj.n, obj.nh);
            obj.filters.Yk_est = nan(obj.ny, 1, obj.nh);

            % Store parameters
            obj.Ts = Ts;
            obj.nu = nu;
            obj.n = n;
            obj.ny = ny;
            obj.T = T;
            obj.nj = nj;
            obj.models = models;
            obj.nh = nh;
            obj.r0 = r0;
            obj.P0 = P0;
            obj.x0 = x0;
            obj.p_seq_g_Yk_init = p_seq_g_Yk_init;
            obj.type = "MKF";
            if label == ""
                label = obj.type;
            end
            obj.label = label;

            if reset
                % Initialize variables
                obj.reset()
            end

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Initialize estimate covariance
            obj.Pkp1 = obj.P0;

            % Initialize state and output estimates
            % Note: At initialization at time k = 0, xkp1_est and
            % ykp1_est represent prior estimates of the states,
            % i.e. x_est(k|k-1) and y_est(k|k-1).
            obj.xkp1_est = obj.x0;
            obj.Pkp1 = obj.P0;
            obj.rk = obj.r0;
            obj.rkm1 = [];

            % Initial values of prior conditional probabilities at k = -1 
            obj.p_seq_g_Yk = obj.p_seq_g_Yk_init;

            % Empty vectors to store values for filter calculations
            % p(y(k)|R(k),Y(k-1))
            obj.p_yk_g_seq_Ykm1 = nan(obj.nh, 1);
            % Pr(r(k)|Y(k-1))
            obj.p_rk_g_Ykm1 = nan(obj.nh, 1);
            % Pr(R(k))
            obj.p_rk_g_rkm1 = nan(obj.nh, 1);
            % Pr(R(k)|Y(k-1))
            obj.p_seq_g_Ykm1 = nan(obj.nh, 1);

            % Reset Kalman filter variables
            obj.filters.Xkp1_est = repmat(obj.xkp1_est, 1, 1, obj.nh);
            obj.filters.Pkp1 = repmat(obj.Pkp1, 1, 1, obj.nh);
            obj.filters.Xk_est = nan(obj.n, 1, obj.nh);
            obj.filters.Pk = nan(obj.n, obj.n, obj.nh);
            obj.filters.Yk_est = nan(obj.ny, 1, obj.nh);

            % At initialization at time k = 0, x_est(k|k)
            % and y_est(k|k) have not yet been computed.
            obj.xk_est = nan(obj.n, 1);
            obj.yk_est = nan(obj.ny, 1);
            obj.Pk = nan(obj.n);

        end
        function KF_update(obj, yk, rk)
        % Update Kalman filter estimates using current 
        % measurement, y(k).
            for f = 1:obj.nh
                m = obj.models{rk(f)};
                [obj.filters.Xk_est(:,:,f), obj.filters.Pk(:,:,f), ...
                    obj.filters.Yk_est(:,:,f), obj.filters.Kf(:,:,f), ...
                    obj.filters.Sk(:,:,f)] = ...
                        kalman_update_f( ...
                            m.C, m.R, ...
                            obj.filters.Xkp1_est(:,:,f), ...
                            obj.filters.Pkp1(:,:,f), ...
                            yk ...
                        );
            end
            % TODO: For testing only - remove later
            assert(~any(isnan(obj.filters.Xk_est), 'all'))
            assert(~any(isnan(obj.filters.Pk), 'all'))
        end
        function KF_predict(obj, uk, rk)
        % Calculate Kalman filter predictions of
        % states at next time instant using current
        % input u(k).
            for f = 1:obj.nh
                m = obj.models{rk(f)};
                [obj.filters.Xkp1_est(:,:,f), obj.filters.Pkp1(:,:,f)] = ...
                    kalman_predict_f( ...
                        m.A, m.B, m.Q, ...
                        obj.filters.Xk_est(:,:,f), ...
                        obj.filters.Pk(:,:,f), ...
                        uk ...
                    );
            end
            % TODO: For testing only - remove later
            assert(~any(isnan(obj.filters.Xkp1_est(:,:,f)), 'all'))
            assert(~any(isnan(obj.filters.Pkp1(:,:,f)), 'all'))
        end
        function MKF_prob_update(obj, yk)
        % % Bayesian updates to conditional probability
        % estimates of each hypothesis

            % Compute Pr(r(k)|r(k-1)) based on Markov transition
            % probability matrix
            obj.p_rk_g_rkm1 = prob_rk(obj.rk, obj.T(obj.rkm1, :)');

            % Loop over each hypothesis
            for f = 1:obj.nh

                % Current system mode for this hypothesis
                r = obj.rk(f);

                % Select system model for this mode
                m = obj.models{r};

                % Output prediction y_f_est(k|k-1)
                yk_pred = m.C * obj.filters.Xkp1_est(:,:,f);

                % Covariance of the output prediction error
                cov_yk = m.C * obj.filters.Pkp1(:,:,f) * m.C' + m.R;

                % Make sure covariance matrix is symmetric
                if ~isscalar(cov_yk)
                    cov_yk = triu(cov_yk.',1) + tril(cov_yk);
                end

                % Probability density of output prediction (assuming a
                % multivariate normal distribution) based on prior
                % estimates computed in previous timestep
                obj.p_yk_g_seq_Ykm1(f) = mvnpdf(yk, yk_pred, cov_yk);

            end
            assert(~any(isnan(obj.p_yk_g_seq_Ykm1)))
            assert(~all(obj.p_yk_g_seq_Ykm1 == 0))

            % Compute Pr(r(k)|Y(k-1)) in current timestep from
            % previous estimate (Pr(r(k-1)|Y(k-1))) and transition
            % probabilities
            obj.p_seq_g_Ykm1 = obj.p_rk_g_rkm1 .* obj.p_seq_g_Yk;

            % Bayesian update of Pr(r(k)|Y(k))
            cond_pds = obj.p_yk_g_seq_Ykm1 .* obj.p_seq_g_Ykm1;
            obj.p_seq_g_Yk = cond_pds ./ sum(cond_pds);
            % Note: above calculation normalizes p_seq_g_Yk so that
            % assert(abs(sum(obj.p_seq_g_Yk) - 1) < 1e-15) % is always true

        end
        function update(obj, yk, uk, rk, rkm1)
        % obj.update(yk, uk, rk)
        % updates the estimates of the multi-model Kalman filter
        % and calculates the predictions of the states and output
        % at the next sample time.
        %
        % Arguments:
        %   obs : struct containing the multi-model Kalman filter
        %       variables (see function mkf_filter).
        %   uk : vector (nu, 1) of system inputs at the current 
        %       sample time.
        %   yk : vector (ny, 1) of system output measurements
        %       at the current sample time.
        %   rk : vector, size (nj, 1)
        %       System modes at current time k.
        %   rkm1 : vector, size (nj, 1) (optional)
        %       System modes at time k - 1. If not specified,
        %       rkm1 is set to the values stored in obj.rk
        %       from the last call to this function.
        %

            % Check size of arguments passed
            assert(isequal(size(uk), [obj.nu 1]), "ValueError: size(uk)")
            assert(isequal(size(yk), [obj.ny 1]), "ValueError: size(yk)")
            assert(isequal(size(rk), [obj.nh 1]), "ValueError: size(rk)")

            if nargin < 5
                % Default: set r(k-1) to previous values of r(k)
                obj.rkm1 = obj.rk;
            else
                % Use provided value for r(k-1)
                obj.rkm1 = rkm1;
            end

            % Store system modes at current time r(k)
            obj.rk = rk;

            % Kalman filter update step
            obj.KF_update(yk, rk)

            % Bayesian updating of hypothesis probabilities
            obj.MKF_prob_update(yk)

            % Merge the Kalman filter estimates and error covariances
            % using a weighted-average based on the hypothesis
            % probabilities
            [obj.xk_est, obj.Pk, obj.yk_est, ~] = ...
                merge_estimates( ...
                    obj.filters.Xk_est, ...
                    obj.filters.Pk, ...
                    obj.filters.Yk_est, ...
                    obj.p_seq_g_Yk ...
                );

            % Kalman filter prediction step - computes estimates
            % of states and error covariances at the next time
            % instant
            obj.KF_predict(uk, rk)

            % TODO: Do we still need a merged xkp1 estimate?
            weights = reshape(obj.p_seq_g_Yk, 1, 1, []);
            obj.xkp1_est = sum(weights .* obj.filters.Xkp1_est, 3);
            assert(~any(isnan(obj.xkp1_est)))
            Xkp1_devs = obj.xkp1_est - obj.filters.Xkp1_est;
            obj.Pkp1 = sum(weights .* (obj.filters.Pkp1 + ...
                pagemtimes(Xkp1_devs, pagetranspose(Xkp1_devs))), 3);

        end
    end
end
