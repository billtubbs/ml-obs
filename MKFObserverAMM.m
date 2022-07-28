% Multi-model observer class definition
%
% obs = MKFObserverAMM(models,Ts,P0,label,x0,y0,p_seq_g_Yk_init)
% Class for simulating an autonomous multi-model (AMM)
% Kalman filter for state estimation of system with more than
% one possible mode (multiple models). This is used as the base
% class for the MKFObserverGPB1 - 1st order Generalised Pseudo-
% Bayes algorithm.
% 
% These observers are in 'filtering form', which means they
% produce posterior estimates of the states and outputs at the 
% current time given the data up to the current time:
%
%  x_est(k|k) : estimate of states at time k
%  y_est(k|k) : estimate of outputs at time k
%
% Predictions (prior estimates) of the states and outputs
% at the next time instant given the data at the current
% time are also calculated:
%
%  x_pred(k+1|k) : estimate of states at time k + 1
%  y_pred(k+1|k) : estimate of outputs at time k + 1
%
% The system model is defined as:
%
%   x(k+1) = A(k) x(k) + B(k) u(k) + w(k)
%     y(k) = C(k) x(k) + v(k)
%
% Note: Assumes no direct transmission (D = 0).
%
% Arguments:
%   models : cell array of structs defining a set of 
%       discrete-time linear system models. Each struct
%       must have the system matrices A, B, C, and the
%       covariance matrices Q and R as fields.
%   Ts : Sample period.
%   P0 : Initial covariance of the state estimation errors
%       (same for each filter).
%   label : string name (optional, default "MKF_AMM").
%   x0 : Initial state estimates (optional, default zeros).
%   y0 : Initial output estimates (optional, default zeros).
%   p_seq_g_Yk_init : (optional, default uniform)
%       Initial prior hypothesis probabilities at time k-1.
%       If not specified, default is equal, i.e. uniform,
%       probability assigned to each hypothesis.
%

classdef MKFObserverAMM < AbstractMKFObserver
    properties
        gamma_k double {mustBeInteger, mustBeNonnegative}
        p_seq_g_Yk_init double
        p_seq_g_Ykm1 double
        p_seq_g_Yk double
        p_yk_g_seq_Ykm1 double
        p_gammak_g_Ykm1 double
        p_gamma_k_g_gamma_km1 double
        Skf (:, :, :) double
        Kf (:, :, :) double
        Pkf (:, :, :) double
        xkp1_est (:, 1) double
        Pkp1 double
        ykp1_est (:, 1) double
        Xkf_est (:, 1, :) double
        Ykf_est (:, 1, :) double
        Pkp1f (:, :, :) double
        Xkp1f_est (:, 1, :) double
        Ykp1f_est (:, 1, :) double
    end
    methods
        function obj = MKFObserverAMM(models,Ts,P0,label,x0,y0, ...
                p_seq_g_Yk_init)

            type = "MKF_AMM";

            % Number of switching systems
            nj = numel(models);

            % Number of states (all models must have same structure)
            n = size(models{1}.A, 1);

            % Number of filters required (= no. of models in
            % the case of AMM observer)
            n_filt = nj;

            if nargin < 7
                % Initial values of prior conditional probabilities at 
                % k = -1. In absence of prior knowledge, assume all 
                % equally likely
                p_seq_g_Yk_init = ones(n_filt, 1) ./ double(n_filt);
            end
            if nargin < 5
                x0 = zeros(n,1);
            end
            if nargin < 6
                y0 = models{1}.C * x0;
            end
            if nargin < 4
                label = type;
            end

            % Call super-class constructor
            obj = obj@AbstractMKFObserver(models,Ts,n_filt,label,x0,y0);

            % Initial value of state estimate error covariance
            obj.P0 = P0;

            % Model indicator values gamma(k) are static - one 
            % filter for each model (zero-based index)
            obj.gamma_k = (0:nj-1)';
            % Transition probabilities Pr(Gamma(k)|Gamma(k-1)) 
            %  = 1 since no transitions with AMM
            obj.p_gamma_k_g_gamma_km1 = ones(obj.nj, 1);

            % Save properties
            obj.p_seq_g_Yk_init = p_seq_g_Yk_init;
            obj.type = type;

            % Initialize variables
            obj.reset();

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Call super-class reset method
            reset@AbstractMKFObserver(obj)

            % Prior conditional probabilities at time k = -1 
            obj.p_seq_g_Yk = obj.p_seq_g_Yk_init;

            % Empty vectors to store variables
            % p(y(k)|Gamma(k),Y(k-1))
            obj.p_yk_g_seq_Ykm1 = nan(obj.n_filt, 1);
            % Pr(gamma(k)|Y(k-1))
            obj.p_gammak_g_Ykm1 = nan(obj.n_filt, 1);
            % Pr(Gamma(k)|Y(k-1))
            obj.p_seq_g_Ykm1 = nan(obj.n_filt, 1);

            % Initialize state and output estimates
            % Note: At initialization at time k = 0, xkp1_est and
            % ykp1_est represent prior estimates of the states,
            % i.e. x_est(k|k-1) and y_est(k|k-1).
            obj.Pkp1 = obj.P0;
            obj.xkp1_est = obj.x0;
            obj.ykp1_est = obj.y0;
            obj.Pkp1f = repmat(obj.P0, 1, 1, obj.n_filt);
            obj.Xkp1f_est = repmat(obj.x0, 1, 1, obj.n_filt);
            obj.Ykp1f_est = repmat(obj.y0, 1, 1, obj.n_filt);

            % Arrays to store variables for each filter
            obj.Skf = nan(obj.ny, obj.ny, obj.n_filt);
            obj.Kf = nan(obj.n, obj.ny, obj.n_filt);
            obj.Pkf = nan(obj.n, obj.n, obj.n_filt);
            obj.Xkf_est = nan(obj.n, 1, obj.n_filt);
            obj.Ykf_est = nan(obj.ny, 1, obj.n_filt);

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk)
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
        %

            assert(isequal(size(uk), [obj.nu 1]), "ValueError: size(uk)")
            assert(isequal(size(yk), [obj.ny 1]), "ValueError: size(yk)")

            % Update all filters
            for f = 1:obj.n_filt

                % Model
                m = obj.models{obj.gamma_k(f) + 1};

                % Update filter estimates and covariance matrix
                [obj.Kf(:,:,f), obj.Xkf_est(:,:,f), obj.Pkf(:,:,f), ...
                    obj.Ykf_est(:,:,f), obj.Skf(:,:,f)] = ...
                    kalman_update_f(m.C, m.R, obj.Xkp1f_est(:,:,f), ...
                        obj.Pkp1f(:,:,f), yk);

                % Prediction step for all filters
                % Kalman filter prediction step
                [obj.Xkp1f_est(:,:,f), obj.Ykp1f_est(:,:,f), ...
                    obj.Pkp1f(:,:,f)] = kalman_predict_f(m.A, m.B, ...
                        m.C, m.Q, obj.Xkf_est(:,:,f), obj.Pkf(:,:,f), uk);

                % Compute posterior probability density of y(k)
                % using the posterior PDF (assumed to be a normal 
                % distribution) and output estimate.

                % Get updated y_f_est(k/k) estimate
                ykf_est = obj.Ykf_est(:,:,f);

                % Make sure covariance matrix is symmetric
                if ~isscalar(obj.Skf(:,:,f))
                    Sk = triu(obj.Skf(:,:,f).',1) + tril(obj.Skf(:,:,f));
                else
                    Sk = obj.Skf(:,:,f);
                end

                % Calculate normal probability density (multivariate)
                obj.p_yk_g_seq_Ykm1(f) = mvnpdf(yk, ykf_est, Sk);

            end

            assert(~any(isnan(obj.p_yk_g_seq_Ykm1)))
            assert(~all(obj.p_yk_g_seq_Ykm1 == 0))

            % Bayesian update to conditional probabilities
            % Compute Pr(Gamma(k)|Y(k-1)) in current timestep from
            % previous estimate (Pr(Gamma(k-1)|Y(k-1))) and transition
            % probabilities
            obj.p_seq_g_Ykm1 = obj.p_gamma_k_g_gamma_km1 .* obj.p_seq_g_Yk;

            % Bayesian update of Pr(Gamma(k)|Y(k))
            likelihood = obj.p_yk_g_seq_Ykm1 .* obj.p_seq_g_Ykm1;
            obj.p_seq_g_Yk = likelihood ./ sum(likelihood);
            % To prevent likelihoods going to zero use this:
            %obj.p_seq_g_Yk = likelihood * 0.998 ./ sum(likelihood) + 0.001;
            % Note: above calculation normalizes p_seq_g_Yk so that
            % assert(abs(sum(obj.p_seq_g_Yk) - 1) < 1e-15) % is always true

            % Compute multi-model observer state and output estimates
            % and estimated state error covariance using the weighted-
            % averages based on the conditional probabilities.
            [obj.xk_est, obj.yk_est, obj.Pk] = ...
                weighted_avg_estimates(obj.Xkf_est, obj.Ykf_est, ...
                    obj.Pkf, obj.p_seq_g_Yk);

            % Weighted average estimates at next time instant
            [obj.xkp1_est, obj.ykp1_est, obj.Pkp1] = ...
                weighted_avg_estimates(obj.Xkp1f_est, obj.Ykp1f_est, ...
                    obj.Pkp1f, obj.p_seq_g_Yk);

        end
    end
end
