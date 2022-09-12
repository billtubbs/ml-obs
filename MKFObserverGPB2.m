% Multi-model Kalman Filter class definition
%
% obs = MKFObserverGPB2(models,P0,T,label,x0,p_seq_g_Yk_init)
% Class for simulating the generalised pseudo-Bayes multi-
% model Kalman filter for state estimation of Markov jump
% linear systems. This is the second-order version of the
% algorithm (GPB2).
%
% This is the filtering form of the observer which 
% produces posterior estimates of the states and
% outputs at the current time instant given the data
% at the current time:
%
%   x_hat(k|k) : estimate of states at time k
%   y_hat(k|k) : estimate of outputs at time k
%
% Prior estimates of the states at the next time
% instant given the data at the current time are also
% calculated:
%
%   x_hat(k+1|k) : estimate of states at time k + 1
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
%   p_seq_g_Yk_init : (optional, default uniform)
%       Initial prior hypothesis probabilities at time k-1.
%       If not specified, default is equal, i.e. uniform,
%       probability assigned to each hypothesis.
%

classdef MKFObserverGPB2 < MKFObserver
    properties
        merged struct
        rkm1 (:, 1) double {mustBeInteger}
    end
    methods
        function obj = MKFObserverGPB2(models,P0,T,label,x0, ...
                p_seq_g_Yk_init)
            arguments
                models (1, :) cell
                P0 double
                T double
                label (1, 1) string = ""
                x0 = []
                p_seq_g_Yk_init = []
            end

            % Get number of system models and check their dimensions
            nj = check_models(models);

            % System modes to be modelled
            r0 = reshape(repmat(1:nj, nj, 1), [], 1);

            % Split prior probabilities into nj^2 copies
            p_seq_g_Yk_init = reshape( ...
                repmat(p_seq_g_Yk_init ./ nj, nj, 1), [], 1);

            % Create super-class observer instance
            obj = obj@MKFObserver(models,P0,T,r0,label,x0,p_seq_g_Yk_init);

            % Store parameters
            obj.type = "MKF_GPB2";

            % Initialize all variables
            obj.reset()

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Call reset method of super class
            reset@MKFObserver(obj);

            % Initialize estimate covariance
            obj.Pkp1 = obj.P0;

            % Create struct to store merged estimates
            obj.merged = struct();
            obj.merged.Xk_est = nan(obj.n, 1, obj.nh);
            obj.merged.Pk = nan(obj.n, obj.n, obj.nh);
            obj.merged.Yk_est = nan(obj.ny, 1, obj.nh);
            obj.merged.p_seq_g_Yk = nan;

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk, rk)
        % updates the estimates of the switching Kalman filter
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
        %

            % Check size of arguments passed
            assert(isequal(size(uk), [obj.nu 1]), "ValueError: size(uk)")
            assert(isequal(size(yk), [obj.ny 1]), "ValueError: size(yk)")

            % Vectors of transitions modelled
            obj.rkm1 = [1 2 1 2]';
            obj.rk = [1 1 2 2]';

            % Update state and output estimates based on current
            % measurement and prior predictions
            [obj.xk_est, obj.yk_est, obj.Pk, obj.merged.Xk_est, ...
                obj.merged.Yk_est, obj.merged.Pk, obj.merged.p_seq_g_Yk] = ...
                GPB2_update( ...
                    obj.models, ...
                    obj.T, ...
                    obj.filters.Xkp1_est, ...
                    obj.filters.Pkp1, ...
                    yk, ...
                    obj.p_seq_g_Yk ...
                );
            assert(~any(isnan(obj.xk_est)))
            assert(~any(isnan(obj.merged.Xk_est)))
            assert(~any(isnan(obj.p_seq_g_Yk)))

            % Calculate predictions of each filter
            % GBP2 branches the estimates from previous time
            % instant when making predictions for next:
            %   xi_est(k+1|k) = Ai(k) * x_est(k|k-1) + Bi(k) * u(k);
            %   Pi(k+1|k) = Ai(k) * P(k|k-1) * Ai(k)' + Qi(k);
            %
            for f = 1:obj.nh
                m = obj.models{obj.rk(f)};
                [obj.filters.Xkp1_est(:,:,f), obj.filters.Pkp1(:,:,f)] = ...
                    kalman_predict_f( ...
                        m.A, m.B, m.Q, ...
                        obj.xk_est, ...
                        obj.Pk, ...
                        uk ...
                    );
            end

            % TODO: Do we need a merged xkp1 estimate?
            weights = reshape(obj.p_seq_g_Yk, 1, 1, []);
            obj.xkp1_est = sum(weights .* obj.filters.Xkp1_est, 3);
            assert(~any(isnan(obj.xkp1_est)))
            Xkp1_devs = obj.xkp1_est - obj.filters.Xkp1_est;
            obj.Pkp1 = sum(weights .* (obj.filters.Pkp1 + ...
                pagemtimes(Xkp1_devs, pagetranspose(Xkp1_devs))), 3);

        end
    end
end