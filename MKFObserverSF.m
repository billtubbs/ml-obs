% Multi-model Kalman Filter class definition
%
% obs = MKFObserverSF(models,P0,seq,T,label,x0,gamma_init, ...
%     p_seq_g_Yk_init)
% Class for simulating a sub-optimal multi-model observer with 
% sequence fusion, for state estimation of a Markov jump linear 
% system. This version produces posterior estimates of the 
% states and outputs at the current time instant given the data
% available at the current time:
%
%  x_hat(k|k) : estimate of states at time k
%  y_hat(k|k) : estimate of outputs at time k
%
% Prior estimates of the states and outputs at the next
% time instant are also calculated:
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
%   seq : model indicator sequences for each filter (in rows).
%   T : Transition probabity matrix for the Markov switching
%       process. T(i,j) represents the probability of the
%       system switching from model i to j.
%   label : string
%       Arbitrary name to identify the observer.
%   x0 : (n, 1) double
%       Initial state estimates (optional, default zeros).
%   gamma_init : (nh, 1) integer (optional, default zeros)
%       Initial prior model indicator value at time k-1 
%       (zero-based, i.e. 0 is for first model).
%   p_seq_g_Yk_init : (nh, 1) double (optional, default uniform)
%       Initial prior probabilities of each hypothesis at
%       time k-1. If not specified, default is equal, i.e.
%       uniform, probability assigned to each hypothesis.
%

classdef MKFObserverSF < matlab.mixin.Copyable
    properties (SetAccess = immutable)
        nu (1, 1) double {mustBeInteger, mustBeNonnegative}
        ny (1, 1) double {mustBeInteger, mustBeNonnegative}
        n (1, 1) double {mustBeInteger, mustBeNonnegative}
        nj (1, 1) double {mustBeInteger, mustBeNonnegative}
        f (1, 1) double {mustBeInteger, mustBeNonnegative}
        nh (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        models cell
        seq (:, 1) cell
        T double
        label (1, 1) string
        P0 double
        x0 (:, 1) double
        gamma_init double {mustBeInteger, mustBeNonnegative}
        p_seq_g_Yk_init double
        i (1, 1) {mustBeInteger, mustBeNonnegative}
        i_next (1, 1) {mustBeInteger, mustBeNonnegative}
        gamma_k double
        p_seq_g_Ykm1 double
        p_seq_g_Yk double
        p_yk_g_seq_Ykm1 double
        p_gammak_g_Ykm1 double
        p_gamma_k double
        filters struct
        idx_branch (1, :) cell
        idx_modes (1, :) cell
        idx_merge (1, :) cell
        xk_est (:, 1) double
        Pk double
        yk_est (:, 1) double
        xkp1_est (:, 1) double
        Pkp1 double
        ykp1_est (:, 1) double
        type (1, 1) string
    end
    methods
        function obj = MKFObserverSF(models,P0,seq,T,label,x0, ...
                gamma_init,p_seq_g_Yk_init)

            % System dimensions
            [n, nu, ny] = check_dimensions(models{1}.A, models{1}.B, ...
                models{1}.C);

            % Number of switching systems
            nj = numel(models);

            % Number of hypotheses
            nh = size(seq, 1);

            if nargin < 8
                % Initial values of prior conditional probabilities at 
                % k = -1. In absence of prior knowledge, assume all 
                % equally likely
                p_seq_g_Yk_init = ones(nh, 1) ./ double(nh);
            end
            if nargin < 7
                % Default assumption about model indicator values at k = -1
                gamma_init = 0;
            end
            if nargin < 6
                x0 = zeros(n,1);
            end

            % Generate the vectors of indices which define the
            % hypothesis branching, mode transitions, and merging 
            % steps of the sequence fusion algorithm
            [obj.idx_branch, obj.idx_modes, obj.idx_merge] = ...
                seq_fusion_indices(cell2mat(seq), nj);

            obj.models = models;
            obj.seq = seq;
            obj.T = T;
            obj.label = label;
            obj.P0 = P0;
            obj.x0 = x0;

            % Fusion horizon length
            obj.f = size(cell2mat(obj.seq), 2);

            % Prior assumptions at initialization
            if isscalar(gamma_init)
                % In case single value specified
                gamma_init = gamma_init * ones(nh, 1);
            end
            obj.gamma_init = gamma_init;
            obj.p_seq_g_Yk_init = p_seq_g_Yk_init;

            % Check transition probability matrix
            % TODO: Is this necessary? If all possible hypotheses 
            % are not modelled (as is the case with some observers)
            % then perhaps T should not be whole?
            % assert(all(abs(sum(obj.T, 2) - 1) < 1e-15), "ValueError: T")
            % TODO: Could calculate beta parameter here, i.e. total 
            % probability measured?


            % Store other useful variables
            obj.nu = nu;
            obj.ny = ny;
            obj.n = n;
            obj.nj = nj;
            obj.nh = nh;
            obj.type = "MKF_SF";

            % Initialize all variables
            obj.reset()

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Switching variable at previous time instant
            obj.gamma_k = obj.gamma_init;

            % Reset sequence index
            obj.i = int16(0);
            obj.i_next = int16(1);

            % Initial values of prior conditional probabilities at 
            % k = -1 
            obj.p_seq_g_Yk = obj.p_seq_g_Yk_init;

            % Empty vectors to store values for filter calculations
            % p(y(k)|Gamma(k),Y(k-1))
            obj.p_yk_g_seq_Ykm1 = nan(obj.nh, 1);
            % Pr(gamma(k)|Y(k-1))
            obj.p_gammak_g_Ykm1 = nan(obj.nh, 1);
            % Pr(Gamma(k))
            obj.p_gamma_k = nan(obj.nh, 1);
            % Pr(Gamma(k)|Y(k-1))
            obj.p_seq_g_Ykm1 = nan(obj.nh, 1);

            % Initialize state and output estimates
            % Note: At initialization at time k = 0, xkp1_est and
            % ykp1_est represent prior estimates of the states,
            % i.e. x_est(k|k-1) and y_est(k|k-1).
            obj.xkp1_est = obj.x0;
            obj.ykp1_est = obj.models{1}.C * obj.xkp1_est;
            obj.Pkp1 = obj.P0;

            % Create struct to store Kalman filter variables
            obj.filters = struct();
            obj.filters.Xkp1_est = repmat(obj.xkp1_est, 1, 1, obj.nh);
            obj.filters.Pkp1 = repmat(obj.P0, 1, 1, obj.nh);
            obj.filters.Ykp1_est = repmat(obj.ykp1_est, 1, 1, obj.nh);

            % At initialization at time k = 0, x_est(k|k)
            % and y_est(k|k) have not yet been computed.
            obj.xk_est = nan(obj.n, 1);
            obj.yk_est = nan(obj.ny, 1);
            obj.Pk = nan(obj.n);

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

            % Increment sequence index
            obj.i = obj.i_next;
            obj.i_next = mod(obj.i, obj.f) + 1;

            % Update model indicator values gamma(k) with the
            % current values from the filter's sequence and keep a
            % copy of the previous values
            gamma_km1 = obj.gamma_k;
            obj.gamma_k = cellfun(@(x) x(:, obj.i), obj.seq);

            % Determine Pr(gamma(k)|gamma(k-1)) based on Markov transition
            % probability matrix
            % TODO: This doesn't need to be a property since gamma_k
            % and p_gamma are properties.
            obj.p_gamma_k = prob_gamma(obj.gamma_k, obj.T(gamma_km1+1, :)');

            % Arrays to collect estimates from each filter
            Xkf_est = zeros(obj.n, 1, obj.nh);
            Pkf_est = zeros(obj.n, obj.n, obj.nh);
            Ykf_est = zeros(obj.ny, 1, obj.nh);
            Xkp1f_est = zeros(obj.n, 1, obj.nh);
            Ykp1f_est = zeros(obj.ny, 1, obj.nh);

            % Bayesian update to conditional probabilities
            for f = 1:obj.nh

                % Compute posterior probability density of y(k)
                % using posterior PDF (normal distribution) and
                % estimates computed in previous timestep

                % Get y_f_est(k/k-1) estimated in previous time step
                yk_est = obj.filters{f}.ykp1_est;

                % Calculate covariance of the output estimation errors
                C = obj.filters{f}.C;
                yk_cov = C * obj.filters{f}.Pkp1 * C' + obj.filters{f}.R;

                % Make sure covariance matrix is symmetric
                if ~isscalar(yk_cov)
                    yk_cov = triu(yk_cov.',1) + tril(yk_cov);
                end

                % Calculate normal probability density (multivariate)
                obj.p_yk_g_seq_Ykm1(f) = mvnpdf(yk, yk_est, yk_cov);

                % Model index at current sample time
                ind = obj.gamma_k(f) + 1;  % MATLAB indexing

                % Set filter model according to current index value
                obj.filters{f}.A = obj.A{ind};
                obj.filters{f}.B = obj.B{ind};
                obj.filters{f}.C = obj.C{ind};
                obj.filters{f}.Q = obj.Q{ind};
                obj.filters{f}.R = obj.R{ind};

                % Update observer estimates, gain and covariance matrix
                obj.filters{f}.update(yk, uk);
                assert(~any(isnan(obj.filters{f}.xkp1_est)))

                % Save state and output estimates
                Xkf_est(:, :, f) = obj.filters{f}.xk_est';
                Pkf_est(:, :, f) = obj.filters{f}.Pk;
                Ykf_est(:, :, f) = obj.filters{f}.yk_est';
                Xkp1f_est(:, :, f) = obj.filters{f}.xkp1_est';
                Pkp1f_est(:, :, f) = obj.filters{f}.Pk;
                Ykp1f_est(:, :, f) = obj.filters{f}.ykp1_est';

            end

            assert(~any(isnan(obj.p_yk_g_seq_Ykm1)))
            assert(~all(obj.p_yk_g_seq_Ykm1 == 0))

            % Compute Pr(Gamma(k)|Y(k-1)) in current timestep from
            % previous estimate (Pr(Gamma(k-1)|Y(k-1))) and transition
            % probabilities
            obj.p_seq_g_Ykm1 = obj.p_gamma_k .* obj.p_seq_g_Yk;

            % Bayesian update of Pr(Gamma(k)|Y(k))
            likelihood = obj.p_yk_g_seq_Ykm1 .* obj.p_seq_g_Ykm1;
            obj.p_seq_g_Yk = likelihood ./ sum(likelihood);
            % Note: above calculation normalizes p_seq_g_Yk so that
            % assert(abs(sum(obj.p_seq_g_Yk) - 1) < 1e-15) % is always true

            % Compute multi-model observer state and output estimates
            % and estimated state error covariance using the weighted-
            % averages based on the conditional probabilities.
            weights = reshape(obj.p_seq_g_Yk, 1, 1, []);
            obj.xk_est = sum(weights .* Xkf_est, 3);
            obj.yk_est = sum(weights .* Ykf_est, 3);
            Xkf_devs = obj.xk_est - Xkf_est;
            obj.Pk = sum(weights .* (Pkf_est + ...
                pagemtimes(Xkf_devs, pagetranspose(Xkf_devs))), 3);
            % TODO: is this correct, how to calculate Pkp1?
            Xkf_devs = obj.xkp1_est - Xkp1f_est;
            obj.Pkp1 = sum(weights .* (Pkp1f_est + ...
                pagemtimes(Xkf_devs, pagetranspose(Xkf_devs))), 3);
            assert(~any(isnan(obj.xk_est)))
            assert(~any(isnan(obj.yk_est)))

            % Weighted averages
            obj.xk_est = sum(weights .* Xkf_est, 3);
            obj.yk_est = sum(weights .* Ykf_est, 3);
            obj.xkp1_est = sum(weights .* Xkp1f_est, 3);
            obj.ykp1_est = sum(weights .* Ykp1f_est, 3);

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
