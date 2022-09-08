% Kalman Filter class definition
%
% obs = KalmanFilterF(A,B,C,Ts,P0,Q,R,label,x0)
% Class for simulating a dynamic Kalman filter
% (i.e. with time-varying gain and estimation error
% covariance). This is the filtering form of the
% KF, which produces posterior estimates of the states and
% outputs at the current time instant given the data
% at the current time:
%
%   x_hat(k|k) : estimate of states at time k
%   y_hat(k|k) : estimate of outputs at time k
%
% Prior estimates of the states and outputs at the next
% time instant given the data at the current time are
% also calculated:
%
%   x_hat(k+1|k) : estimate of states at time k + 1
%   y_hat(k+1|k) : estimate of outputs at time k + 1
%
% The system model is defined as:
%
%   x(k+1) = A(k) x(k) + B(k) u(k) + w(k)
%     y(k) = C(k) x(k) + v(k)
%
% Note: there is no direct transmission (D = 0).
%
% Arguments:
%   A, B, C : matrices
%       Discrete-time system model matrices.
%   Ts : double
%       Sampling period.
%   P0 : matrix, size (n, n)
%       Initial value of covariance matrix of the state
%       estimates at time k = 0.
%   Q : matrix, size (n, n)
%       State error covariance.
%   R : matrix, size (ny, ny)
%       Output measurement noise covariance.
%   label : string (optional)
%       Name.
%   x0 : vector, size(n, 1), (optional)
%       Initial state estimates at time k = 0.
%

classdef KalmanFilterF < matlab.mixin.Copyable
    properties (SetAccess = immutable)
        n (1, 1) double {mustBeInteger, mustBeNonnegative}
        nu (1, 1) double {mustBeInteger, mustBeNonnegative}
        ny (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        model struct
        Ts (1, 1) double {mustBeNonnegative}
        P0 double
        xk_est (:, 1) double
        yk_est (:, 1) double
        xkp1_est (:, 1) double
        ykp1_est (:, 1) double
        Kf double
        Pk double
        Pkp1 double
        Sk double
        label (1, 1) string
        x0 (:, 1) double
        type (1, 1) string  % TODO: is this still needed? use classdef
    end
    methods
        function obj = KalmanFilterF(model,P0,label,x0)
            arguments
                model struct
                P0 double
                label (1, 1) string = ""
                x0 (:, 1) double = []
            end

            % Check system model dimensions
            [n, nu, ny] = check_dimensions(model.A, model.B, model.C);

            % Check size of other parameters
            assert(isequal(size(model.Q), [n n]), "ValueError: size(model.Q)")
            assert(isequal(size(model.R), [ny ny]), "ValueError: size(model.R)")
            assert(isequal(size(P0), [n n]), "ValueError: size(P0)")

            % Determine initial state values
            if isempty(x0)
                x0 = zeros(n, 1);  % default initial states
            else
                assert(isequal(size(x0), [n 1]), "ValueError: size(x0)")
            end

            % Store parameters
            obj.Ts = model.Ts;
            obj.nu = nu;
            obj.n = n;
            obj.ny = ny;
            obj.model = model;
            obj.P0 = P0;
            obj.x0 = x0;
            obj.type = "KFF";
            if label == ""
                label = obj.type;
            end
            obj.label = label;

            % Initialize variables
            obj.reset()

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Initialize state and output estimates
            % Note: At initialization time k = 0, xkp1_est and
            % ykp1_est represent prior estimates of the states,
            % i.e. x_est(k|k-1) and y_est(k|k-1).
            obj.xkp1_est = obj.x0; 
            obj.ykp1_est = obj.model.C * obj.xkp1_est;

            % Initialize error covariance P(k|k-1)
            obj.Pkp1 = obj.P0;

            % At initialization time k = 0, x_est(k|k)
            % and y_est(k|k) have not yet been computed.
            obj.xk_est = nan(obj.n, 1);
            obj.yk_est = nan(obj.ny, 1);
            obj.Pk = nan(obj.n);
            obj.Sk = nan(obj.ny);

            % Gain will be calculated at first update
            obj.Kf = nan(obj.n, obj.ny);

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk) updates the gain and covariance matrix
        % of the Kalman filter and calculates the estimates of the
        % states and output at the current sample time given the 
        % current measurement and inputs.
        %
        % After calling this method, obs.xk_est and obs.yk_est
        % will be the 'a posteriori' estimates of x(k|k) and
        % y(k|k) i.e. based on the data up to time k.
        %
        % Arguments:
        %   yk : vector, size (ny, 1)
        %       System output measurements at current time k.
        %   uk : vector, size (nu, 1), optional
        %       System inputs at current time k.
        %

            % Update estimates based on current measurement
            [obj.xk_est, obj.Pk, obj.yk_est, obj.Kf, obj.Sk] = ...
                kalman_update_f(obj.model.C, obj.model.R, ...
                obj.xkp1_est, obj.Pkp1, yk);

            % Predict states at next time instant
            [obj.xkp1_est, obj.Pkp1] = kalman_predict_f(obj.model.A, ...
                obj.model.B, obj.model.Q, obj.xk_est, obj.Pk, uk);

            % Predicted output at next time instant
            obj.ykp1_est = obj.model.C * obj.xkp1_est;

        end
    end
end