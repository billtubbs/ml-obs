% Multi-model Kalman Filter class definition
%
% obs = MKFObserverF(A,B,C,D,Ts,P0,Q,R,seq,T,label,x0,gamma_init, ...
%     p_seq_g_Yk_init)
% Class for simulating a multi-model Kalman filter for state
% estimation of a Markov jump linear system. This is the
% filtering form of the MKFObserver, which produces posterior
% estimates of the states and outputs at the current time 
% instant given the data at the current time:
%
%  x_hat(k|k) : estimate of states at time k
%  y_hat(k|k) : estimate of outputs at time k
%
% Arguments:
%   A, B, C, D : Cell arrays containing discrete-time system
%       matrices for each switching system modelled.
%   Ts : Sample period.
%   P0 : Initial covariance matrix of the state estimates
%       (same for each filter).
%   Q : Cell array of process noise covariance matrices for
%       each switching system.
%   R : Cell array of output measurement noise covariance
%       matrices for each switching system.
%   seq : model indicator sequences for each filter (in rows).
%   T : Transition probabity matrix of the Markov switching
%       process.
%   label : string name.
%   x0 : Initial state estimates (optional, default zeros).
%   gamma_init : (optional, default zeros)
%       Initial prior model indicator value at time k-1 
%       (zero-based, i.e. 0 is for first model).
%   p_seq_g_Yk_init : (optional, default uniform)
%       Initial prior hypothesis probabilities at time k-1.
%       If not specified, default is equal, i.e. uniform,
%       probability assigned to each hypothesis.
%

classdef MKFObserverF < matlab.mixin.Copyable
    properties (SetAccess = immutable)
        Ts (1, 1) double {mustBeNonnegative}
        n (1, 1) double {mustBeInteger, mustBeNonnegative}
        nu (1, 1) double {mustBeInteger, mustBeNonnegative}
        ny (1, 1) double {mustBeInteger, mustBeNonnegative}
        nj (1, 1) double {mustBeInteger, mustBeNonnegative}
        f (1, 1) double {mustBeInteger, mustBeNonnegative}
        n_filt (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        A cell
        B cell
        C cell
        D cell
        Q cell
        R cell
        seq cell
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
        filters  cell
        xk_est (:, 1) double
        P double
        yk_est (:, 1) double
        type (1, 1) string
    end
    methods
        function obj = MKFObserver(A,B,C,D,Ts,P0,Q,R,seq,T,label,x0, ...
                gamma_init,p_seq_g_Yk_init)

            % System dimensions
            [n, nu, ny] = check_dimensions(A{1}, B{1}, C{1}, D{1});

            % Number of switching systems
            nj = numel(A);

            % Number of filters required
            n_filt = size(seq, 1);

            if nargin < 14
                % Initial values of prior conditional probabilities at 
                % k = -1. In absence of prior knowledge, assume all 
                % equally likely
                p_seq_g_Yk_init = ones(n_filt, 1) ./ double(n_filt);
            end
            if nargin < 13
                % Default assumption about model indicator values at k = -1
                gamma_init = 0;
            end
            if nargin < 12
                x0 = zeros(n,1);
            end
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.D = D;
            obj.Ts = Ts;
            obj.Q = Q;
            obj.R = R;
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
                gamma_init = gamma_init * ones(n_filt, 1);
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

            % Check all other system matrix dimensions have same 
            % input/output dimensions and number of states.
            for j = 2:nj
                [n_j, nu_j, ny_j] = check_dimensions(A{j}, B{j}, C{j}, D{j});
                assert(isequal([n_j, nu_j, ny_j], [n, nu, ny]), ...
                    "ValueError: size of A, B, C, and D")
            end

            % Store other useful variables
            obj.n = n;
            obj.nu = nu;
            obj.ny = ny;
            obj.nj = nj;
            obj.n_filt = n_filt;
            obj.type = "MKF";

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

            % Initial values of prior conditional probabilities at k = -1 
            obj.p_seq_g_Yk = obj.p_seq_g_Yk_init;

            % Empty vectors to store values for filter calculations
            % p(y(k)|Gamma(k),Y(k-1))
            obj.p_yk_g_seq_Ykm1 = nan(obj.n_filt, 1);
            % Pr(gamma(k)|Y(k-1))
            obj.p_gammak_g_Ykm1 = nan(obj.n_filt, 1);
            % Pr(Gamma(k))
            obj.p_gamma_k = nan(obj.n_filt, 1);
            % Pr(Gamma(k)|Y(k-1))
            obj.p_seq_g_Ykm1 = nan(obj.n_filt, 1);

            % Create multi-model observer
            obj.filters = cell(obj.n_filt, 1);
            fmt = strcat('%s%0', ...
                char(string(strlength(string(obj.n_filt)))), 'd');
            % Initialize each filter
            for i = 1:obj.n_filt
                label_i = sprintf(fmt,obj.label,i);
                % Index of system model
                ind = obj.gamma_k(i) + 1;
                obj.filters{i} = KalmanFilterF(obj.A{ind},obj.B{ind}, ...
                    obj.C{ind},obj.D{ind},obj.Ts,obj.P0, obj.Q{ind}, ...
                    obj.R{ind},label_i,obj.x0);
            end

            % Initialize state and output estimates
            % Note: At initialization at time k = 0, x_est(k|k)
            % and y_est(k|k) have not yet been computed.
            obj.xkp1_est = nan(obj.n, 1);
            obj.ykp1_est = nan(obj.ny, 1);

            % At initialization at time k = 0, xkp1_est and
            % ykp1_est represent prior estimates, i.e.
            % x_est(k|k-1) and y_est(k|k-1).
            obj.xkp1_est = obj.x0;
            obj.ykp1_est = obj.C{1} * obj.xkp1_est;

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

            % Compute Pr(gamma(k)|gamma(k-1)) based on Markov transition
            % probability matrix
            % TODO: This doesn't need to be a property since gamma_k
            % and p_gamma are properties.
            obj.p_gamma_k = prob_gamma(obj.gamma_k, obj.T(gamma_km1+1, :)');

            % Arrays to collect estimates from each filter
            Xkf_est = zeros(obj.n, 1, obj.n_filt);
            Pkf_est = zeros(obj.n, obj.n, obj.n_filt);
            Ykf_est = zeros(obj.ny, 1, obj.n_filt);

            % Bayesian update to conditional probabilities
            for f = 1:obj.n_filt

                % Compute posterior probability density of y(k)
                % using posterior PDF (normal distribution) and
                % estimates computed in previous timestep

                % Update observer estimates, gain and covariance matrix
                obj.filters{f}.update(yk, uk);
                assert(~any(isnan(obj.filters{f}.xk_est)))

                % Save state and output estimates for next timestep
                Xkf_est(:, :, f) = obj.filters{f}.xk_est';
                Pkf_est(:, :, f) = obj.filters{f}.P;
                Ykf_est(:, :, f) = obj.filters{f}.yk_est';

                % Get y_f_est(k/k-1) estimated in previous time step
                yk_est = obj.filters{f}.ykp1_est;

                % Calculate covariance of the output estimation errors
                C = obj.filters{f}.C;
                yk_cov = C * obj.filters{f}.P * C' + obj.filters{f}.R;

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
                obj.filters{f}.D = obj.D{ind};
                obj.filters{f}.Q = obj.Q{ind};
                obj.filters{f}.R = obj.R{ind};

            end

            assert(~any(isnan(obj.p_yk_g_seq_Ykm1)))
            assert(~all(obj.p_yk_g_seq_Ykm1 == 0))

            % Compute Pr(Gamma(k)|Y(k-1)) in current timestep from
            % previous estimate (Pr(Gamma(k-1)|Y(k-1))) and transition
            % probabilities
            obj.p_seq_g_Ykm1 = obj.p_gamma_k .* obj.p_seq_g_Yk;

            % Bayesian update of Pr(Gamma(k)|Y(k))
            cond_pds = obj.p_yk_g_seq_Ykm1 .* obj.p_seq_g_Ykm1;
            obj.p_seq_g_Yk = cond_pds ./ sum(cond_pds);
            % Note: above calculation normalizes p_seq_g_Yk so that
            % assert(abs(sum(obj.p_seq_g_Yk) - 1) < 1e-15) % is always true

            % Compute multi-model observer state and output estimates
            % and estimated state error covariance using the weighted-
            % averages based on the conditional probabilities.
            weights = reshape(obj.p_seq_g_Yk, 1, 1, []);
            obj.xkp1_est = sum(weights .* Xkf_est, 3);
            obj.ykp1_est = sum(weights .* Ykf_est, 3);
            Xkf_devs = obj.xkp1_est - Xkf_est;
            obj.P = sum(weights .* (Pkf_est + ...
                pagemtimes(Xkf_devs, pagetranspose(Xkf_devs))), 3);
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
