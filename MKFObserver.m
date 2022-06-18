% Multi-model Kalman Filter class definition
%
% Class for simulating a multi-model Kalman filter for state
% estimation of a Markov jump linear system.
%

classdef MKFObserver < matlab.mixin.Copyable
    properties (SetAccess = immutable)
        Ts (1, 1) double {mustBeNonnegative}
        n (1, 1) double {mustBeInteger, mustBeNonnegative}
        nu (1, 1) double {mustBeInteger, mustBeNonnegative}
        ny (1, 1) double {mustBeInteger, mustBeNonnegative}
        nj (1, 1) double {mustBeInteger, mustBeNonnegative}
        d (1, 1) double {mustBeInteger, mustBeNonnegative}
        f (1, 1) double {mustBeInteger, mustBeNonnegative}
        n_filt (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        A cell
        B cell
        C cell
        D cell
        P0 double
        Q cell
        R cell
        seq cell
        gamma0 double {mustBeInteger, mustBeNonnegative}
        T double
        label (1, 1) string
        x0 (:, 1) double
        i (1, 2) {mustBeInteger, mustBeNonnegative}
        i_next (1, 2) {mustBeInteger, mustBeNonnegative}
        gamma_k double
        p_seq_g_Yk double
        p_yk_g_seq_Ykm1 double
        p_gammak_g_Ykm1 double
        p_gamma_k double
        p_seq_g_Ykm1 double
        filters  cell
        xkp1_est (:, 1) double
        ykp1_est (:, 1) double
        type (1, 1) string
    end
    methods
        function obj = MKFObserver(A,B,C,D,Ts,P0,Q,R,seq,T,d,label,x0,gamma0)
        % obs = MKFObserver(A,B,C,D,Ts,P0,Q,R,seq,T,d,label,x0,gamma0)
        %
        % Arguments:
        %	A, B, C, D : cell arrays containing discrete-time system
        %       matrices for each switching system modelled.
        %   Ts : sample period.
        %   P0 : Initial covariance matrix of the state estimates
        %       (same for each filter).
        %   Q : cell array of process noise covariance matrices for
        %       each switching system.
        %   R : cell array of output measurement noise covariance
        %       matrices for each switching system.
        %   seq : model indicator sequences for each filter (in rows).
        %   T : transition probabity matrix of the Markov switching
        %       process.
        %   d : detection interval length in number of sample periods.
        %   label : string name.
        %   x0 : intial state estimates (optional, default zeros)
        %   gamma0 : initial model indicator value (zero-based)
        %       (optional, default zeros).

            % Number of switching systems
            nj = numel(A);

            % System dimensions
            [n, nu, ny] = check_dimensions(A{1}, B{1}, C{1}, D{1});
            if nargin < 14
                gamma0 = 0;
            end
            if nargin < 13
                x0 = zeros(n,1);
            end
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.D = D;
            obj.Ts = Ts;
            obj.P0 = P0;
            obj.Q = Q;
            obj.R = R;
            obj.seq = seq;
            obj.T = T;
            obj.d = d;
            obj.label = label;
            obj.x0 = x0;

            % Number of filters required
            obj.n_filt = size(obj.seq, 1);

            % Assumption about initial model indicator
            if isscalar(gamma0)
                % In case single value specified
                gamma0 = gamma0 * ones(obj.n_filt, 1);
            end
            obj.gamma0 = gamma0;

            % Fusion horizon length
            obj.f = size(cell2mat(obj.seq), 2);

            % Check transition probability matrix
            assert(all(abs(sum(obj.T, 2) - 1) < 1e-15), "ValueError: T")

            % Check all other system matrix dimensions have same 
            % input/output dimensions and number of states.
            for j = 2:nj
                [n_j, nu_j, ny_j] = check_dimensions(A{j}, B{j}, C{j}, D{j});
                assert(isequal([n_j, nu_j, ny_j], [n, nu, ny]), ...
                    "ValueError: size of A, B, C, and D")
            end

            % Initialize all variables
            obj.reset()

            % Add other useful variables
            obj.n = n;
            obj.nu = nu;
            obj.ny = ny;
            obj.nj = nj;
            obj.type = "MKF";

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %
            % Switching variable at previous time instant
            obj.gamma_k = obj.gamma0;

            % Sequence index and counter for prob. updates
            % obj.i(1) is the sequence index (1 <= i(1) <= obj.f)
            % obj.i(2) is the counter for prob. updates (1 <= i(2) <= obj.d)
            obj.i = int16([0 0]);
            obj.i_next = int16([1 1]);

            % Initialize conditional probabilities: all equally likely
            obj.p_seq_g_Yk = ones(obj.n_filt, 1) ./ double(obj.n_filt);

            % Empty vectors to store values for filter calculations
            % p(y(k)|Gamma(k),Y(k-1))
            obj.p_yk_g_seq_Ykm1 = zeros(obj.n_filt, 1);
            % Pr(gamma(k)|Y(k-1))
            obj.p_gammak_g_Ykm1 = zeros(obj.n_filt, 1);
            % Pr(gamma(k))
            obj.p_gamma_k = zeros(obj.n_filt, 1);
            % Pr(Gamma(k)|Y(k-1))
            obj.p_seq_g_Ykm1 = zeros(obj.n_filt, 1);

            % Create multi-model observer
            obj.filters = cell(obj.n_filt, 1);
            fmt = strcat('%s%0', ...
                char(string(strlength(string(obj.n_filt)))), 'd');
            % Initialize each filter
            for i = 1:obj.n_filt
                label_i = sprintf(fmt,obj.label,i);
                % Index of system model
                ind = obj.gamma_k(i) + 1;
                obj.filters{i} = KalmanFilter(obj.A{ind},obj.B{ind}, ...
                    obj.C{ind},obj.D{ind},obj.Ts,obj.P0, obj.Q{ind}, ...
                    obj.R{ind},label_i,obj.x0);
            end

            % Initialize estimates
            obj.xkp1_est = obj.x0;
            obj.ykp1_est = obj.C{1} * obj.xkp1_est;

        end
        function update(obj, yk, uk, show_plots)
        % obj.update(yk, uk, show_plots)
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

            % Debugging option
            if nargin < 4
                show_plots = false;
            end

            assert(isequal(size(uk), [obj.nu 1]), "ValueError: size(uk)")
            assert(isequal(size(yk), [obj.ny 1]), "ValueError: size(yk)")

            % Increment sequence index and update counter
            % obj.i(1) is the sequence index (1 <= i(1) <= obj.f)
            % obj.i(2) is the counter for prob. updates (1 <= i(2) <= obj.d)
            % Whenever obj.i(2) exceeds obj.d (the spacing parameter), it is
            % reset to 1, and the sequence index obj.i(1) is incremented.
            obj.i = obj.i_next;
            obj.i_next = [mod(obj.i(1) - 1 + ...
                          idivide(obj.i(2), obj.d), obj.f) + 1, ...
                          mod(obj.i(2), obj.d) + 1];

            % Arrays to store model indicator sequence values
            obj.p_gamma_k = nan(obj.n_filt, 1);

            % Arrays to collect estimates from each filter
            Xkf_est = zeros(obj.n_filt, obj.n);
            Ykf_est = zeros(obj.n_filt, obj.ny);

            % Bayesian update to conditional probabilities
            for f = 1:obj.n_filt

                % Compute posterior probability density of y(k)
                % using posterior PDF (normal distribution) and
                % estimates computed in previous timestep

                % Get y_est(k/k-1) estimated in previous time step
                yk_est = obj.filters{f}.ykp1_est;

                % Calculate covariance of the output estimation errors
                P = obj.filters{f}.P;
                C = obj.filters{f}.C;
                yk_cov = C*P*C' + obj.filters{f}.R;

                % Make sure covariance matrix is symmetric
                if ~isscalar(yk_cov)
                    yk_cov = triu(yk_cov.',1) + tril(yk_cov);
                end

                % Calculate normal probability density (multivariate)
                obj.p_yk_g_seq_Ykm1(f) = mvnpdf(yk, yk_est, yk_cov);

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
                    legend('$p(y(k))$', '$\hat{y}(k)$', '$y_m(k)$', ...
                        'Interpreter','Latex')
                    set(gcf,'Position',[f*250-150 50 250 150])
                end

                % Index of model used in previous sample time
                ind_km1 = obj.gamma_k(f) + 1;  % MATLAB indexing

                % Update model indicator value gamma(k) with the
                % current value from the filter's sequence
                obj.gamma_k(f) = obj.seq{f}(:, obj.i(1));

                % Model index at current sample time
                ind = obj.gamma_k(f) + 1;  % MATLAB indexing

                % Compute Pr(gamma(k)) based on Markov transition
                % probability matrix
                % TODO: Can this be computed prior to the for loop for 
                % efficiency?
                obj.p_gamma_k(f) = prob_gamma(obj.gamma_k(f), ...
                    obj.T(ind_km1, :)');

                % Update filter covariances if at start of a detection
                % interval.
                % Note: Currently, this update must happen every sample 
                % period so that S-functions do not have to memorize all
                % the model parameters each timestep.
                %if obj.i(2) == 1
                % Select filter system model based on current
                % model indicator value
                % TODO: only need to change these if ind changed?
                obj.filters{f}.A = obj.A{ind};
                obj.filters{f}.B = obj.B{ind};
                obj.filters{f}.C = obj.C{ind};
                obj.filters{f}.D = obj.D{ind};
                obj.filters{f}.Q = obj.Q{ind};
                obj.filters{f}.R = obj.R{ind};
                %end

                % Update observer estimates, gain and covariance matrix
                obj.filters{f}.update(yk, uk);
                assert(~any(isnan(obj.filters{f}.xkp1_est)))

                % Save state and output estimates for next timestep
                Xkf_est(f, :) = obj.filters{f}.xkp1_est';
                Ykf_est(f, :) = obj.filters{f}.ykp1_est';

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

            if show_plots
                % Plot gamma_k, p_yk_g_seq_Ykm1 and 
                figure(50)
                subplot(3,1,1)
                bar(1:obj.n_filt, obj.p_yk_g_seq_Ykm1)
                title('$p(y(k)|\Gamma(k),Y(k-1))$','Interpreter','Latex')
                subplot(3,1,2)
                bar(1:obj.n_filt, obj.p_gamma_k)
                ylim([0, 1])
                title('Shock probabilities $\gamma(k)$','Interpreter','Latex')
                subplot(3,1,3)
                bar(1:obj.n_filt, obj.p_seq_g_Yk)
                xlabel('Filter')
                title('$p(\Gamma(k)|Y(k))$','Interpreter','Latex')
            end

            % Compute multi-model observer state and output estimates
            obj.xkp1_est = sum(Xkf_est .* obj.p_seq_g_Yk, 1)';
            obj.ykp1_est = sum(Ykf_est .* obj.p_seq_g_Yk, 1)';
            assert(~any(isnan(obj.xkp1_est)))
            assert(~any(isnan(obj.ykp1_est)))
            % TODO: Calculate multi-model state covariance estimate

            if show_plots 
                disp(obj.label);
                disp('stop')
            end
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
