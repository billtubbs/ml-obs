% Multi-model Kalman Filter class definition
%
% obs = MKFObserver(models,P0,T,r0,label,x0,p_seq_g_Yk_init)
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
%       Initial prior system mode at time k-1 (zero-based,
%       i.e. 0 is for first model).
%   p_seq_g_Yk_init : (optional, default uniform)
%       Initial prior hypothesis probabilities at time k-1.
%       If not specified, default is equal, i.e. uniform,
%       probability assigned to each hypothesis.
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
        K double
        P0 double
        T double
        r0 (:, 1) double {mustBeInteger}
        label (1, 1) string
        x0 (:, 1) double
        p_seq_g_Yk_init double
        p_seq_g_Ykm1 double
        p_seq_g_Yk double
        p_yk_g_seq_Ykm1 double
        p_rk_g_Ykm1 double
        p_rk double
        filters  struct
        xk_est (:, 1) double
        Pk double
        yk_est (:, 1) double
        xkp1_est (:, 1) double
        Pkp1 double
        ykp1_est (:, 1) double
        rk (:, 1) double {mustBeInteger}
        type (1, 1) string
    end
    methods
        function obj = MKFObserver(models,P0,T,r0,label,x0,p_seq_g_Yk_init)
            arguments
                models (1, :) cell
                P0 double
                T double
                r0
                label (1, 1) string = ""
                x0 = []
                p_seq_g_Yk_init = []
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
                p_seq_g_Yk_init = ones(nh, 1) ./ double(nh);
            end

            % Determine initial state values
            if isempty(x0)
                x0 = zeros(obj.n, 1);  % default initial states
            else
                assert(isequal(size(x0), [n 1]))
            end

            % Check transition probability matrix
            assert(isequal(size(T), [nj nj]), "ValueError: size(T)")
            assert(all(abs(sum(obj.T, 2) - 1) < 1e-15), "ValueError: T")

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

            % Initialize all variables
            obj.reset()

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Initialize estimate covariance
            obj.Pkp1 = obj.P0;

            % Gain will be calculated dynamically
            obj.K = nan(obj.n, 1);

            % Initialize state and output estimates
            % Note: At initialization at time k = 0, xkp1_est and
            % ykp1_est represent prior estimates of the states,
            % i.e. x_est(k|k-1) and y_est(k|k-1).
            obj.xkp1_est = obj.x0;
            obj.rk = obj.r0;
            obj.ykp1_est = obj.models{obj.rk(1)}.C * obj.xkp1_est;
            obj.Pkp1 = obj.P0;

            % Initial values of prior conditional probabilities at k = -1 
            obj.p_seq_g_Yk = obj.p_seq_g_Yk_init;

            % Empty vectors to store values for filter calculations
            % p(y(k)|R(k),Y(k-1))
            obj.p_yk_g_seq_Ykm1 = nan(obj.nh, 1);
            % Pr(r(k)|Y(k-1))
            obj.p_rk_g_Ykm1 = nan(obj.nh, 1);
            % Pr(R(k))
            obj.p_rk = nan(obj.nh, 1);
            % Pr(R(k)|Y(k-1))
            obj.p_seq_g_Ykm1 = nan(obj.nh, 1);

            % Create struct to store Kalman filter variables
            obj.filters = struct();
            obj.filters.Xkp1_est = repmat(obj.xkp1_est, 1, 1, obj.nh);
            obj.filters.Pkp1 = repmat(obj.P0, 1, 1, obj.nh);
            obj.filters.Ykp1_est = repmat(obj.ykp1_est, 1, 1, obj.nh);
            obj.filters.Kf = nan(obj.ny, obj.n, obj.nh);

            % At initialization at time k = 0, x_est(k|k)
            % and y_est(k|k) have not yet been computed.
            obj.xk_est = nan(obj.n, 1);
            obj.yk_est = nan(obj.ny, 1);
            obj.Pk = nan(obj.n);

        end
        function update(obj, yk, uk, rk)
        % obj.update(yk, uk, rk)
        % updates the multi-model Kalman filter and calculates the
        % estimates of the states and output at the next sample
        % time.
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
        %

            assert(isequal(size(uk), [obj.nu 1]), "ValueError: size(uk)")
            assert(isequal(size(yk), [obj.ny 1]), "ValueError: size(yk)")

            % Update model indicator values r(k) with the
            % current values from the filter's sequence and keep a
            % copy of the previous values
            rkm1 = obj.rk;
            obj.rk = rk;

            % Compute Pr(r(k)|r(k-1)) based on Markov transition
            % probability matrix
            % TODO: This doesn't need to be a property since rk
            % and p_rk are properties.
            obj.p_rk = prob_rk(obj.rk, obj.T(rkm1, :)');

            % Bayesian update to conditional probabilities
            for f = 1:obj.nh

                % Compute posterior probability density of y(k)
                % using posterior PDF (normal distribution) and
                % estimates computed in previous timestep

                % Get y_f_est(k/k-1) estimated in previous time step
                yk_est = obj.filters.Ykp1_est(:,:,f);

                % System mode
                r = obj.rk(f);

                % Calculate covariance of the output estimation errors
                C = obj.models{r}.C;
                R = obj.models{r}.R;
                yk_cov = C * obj.filters.Pkp1(:,:,f) * C' + R;

                % Make sure covariance matrix is symmetric
                if ~isscalar(yk_cov)
                    yk_cov = triu(yk_cov.',1) + tril(yk_cov);
                end

                % Calculate normal probability density (multivariate)
                obj.p_yk_g_seq_Ykm1(f) = mvnpdf(yk, yk_est, yk_cov);

                % Update correction gain and covariance matrix
                [obj.filters.Kf(:,:,f), obj.filters.Pkp1(:,:,f)] = ...
                    kalman_update(obj.filters.Pkp1(:,:,f), ...
                    obj.models{r}.A, obj.models{r}.C, ...
                    obj.models{r}.Q, obj.models{r}.R);

                % Update state and output estimates in next timestep
                obj.filters.Xkp1_est(:,:,f) = obj.models{r}.A * obj.filters.Xkp1_est(:,:,f) ...
                    + obj.models{r}.B * uk + ...
                    obj.filters.Kf(:,:,f) * (yk - obj.models{r}.C * obj.filters.Xkp1_est(:,:,f));
                % TODO: This is wrong because model in next timestep should
                % be used
                obj.filters.Ykp1_est(:,:,f) = obj.models{r}.C * obj.filters.Xkp1_est(:,:,f);

            end

            assert(~any(isnan(obj.p_yk_g_seq_Ykm1)))
            assert(~all(obj.p_yk_g_seq_Ykm1 == 0))

            % Compute Pr(Gamma(k)|Y(k-1)) in current timestep from
            % previous estimate (Pr(Gamma(k-1)|Y(k-1))) and transition
            % probabilities
            obj.p_seq_g_Ykm1 = obj.p_rk .* obj.p_seq_g_Yk;

            % Bayesian update of Pr(Gamma(k)|Y(k))
            cond_pds = obj.p_yk_g_seq_Ykm1 .* obj.p_seq_g_Ykm1;
            obj.p_seq_g_Yk = cond_pds ./ sum(cond_pds);
            % Note: above calculation normalizes p_seq_g_Yk so that
            % assert(abs(sum(obj.p_seq_g_Yk) - 1) < 1e-15) % is always true

            % Compute multi-model observer state and output estimates
            % and estimated state error covariance using the weighted-
            % averages based on the conditional probabilities.
            weights = reshape(obj.p_seq_g_Yk, 1, 1, []);
            obj.xkp1_est = sum(weights .* obj.filters.Xkp1_est, 3);
            obj.ykp1_est = sum(weights .* obj.filters.Ykp1_est, 3);
            Xkp1_devs = obj.xkp1_est - obj.filters.Xkp1_est;
            obj.Pkp1 = sum(weights .* (obj.filters.Pkp1 + ...
                pagemtimes(Xkp1_devs, pagetranspose(Xkp1_devs))), 3);
            assert(~any(isnan(obj.xkp1_est)))
            assert(~any(isnan(obj.ykp1_est)))

        end
    end
% TODO: Should do a deep copy automatically (i.e. without this)
% - https://www.mathworks.com/help/matlab/matlab_oop/custom-copy-behavior.html
%     methods (Access = protected)
%         function cp = copyElement(obj)
%             % Shallow copy object
%             cp = copyElement@matlab.mixin.Copyable(obj);
%            cp.Prop1 = obj.Prop1;
%            cp.Prop2 = datestr(now);
%         end
%     end
end
