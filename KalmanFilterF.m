% Kalman Filter class definition
%
% obs = KalmanFilterF(A,B,C,D,Ts,P0,Q,R,label,x0)
% Class for simulating a dynamic Kalman filter
% (i.e. with time-varying gain and estimation error
% covariance). This is the filtering form of the
% KF which produces posterior estimates of the states and
% outputs at the current time instant given the data
% at the current time:
%
%  x_hat(k|k) : estimate of states at time k
%  y_hat(k|k) : estimate of outputs at time k
%
% Arguments:
%   A, B, C, D : matrices
%       Discrete-time system model matrices.
%   Ts : double
%       Sampling period.
%   P0 : matrix, size (n, n)
%       Initial value of covariance matrix of the state
%       estimates.
%   Q : matrix, size (n, n)
%       Process noise covariance.
%   R : matrix, size (ny, ny)
%       Output measurement noise covariance.
%   label : string (optional)
%       Name.
%   x0 : vector, size(n, 1), (optional)
%       Initial state estimates.
%

classdef KalmanFilterF < AbstractLinearFilter
    properties
        xk_est (:, 1) double
        yk_est (:, 1) double
        xkp1_est (:, 1) double
        ykp1_est (:, 1) double
        Kf double
        P double
        Q double
        R double
        P0 double
    end
    methods
        function obj = KalmanFilterF(A,B,C,D,Ts,P0,Q,R,varargin)

            % Call super-class constuctor
            obj = obj@AbstractLinearFilter(A,B,C,D,Ts,varargin{:});
            n = obj.n;
            ny = obj.ny;

            % Set additional properties for dynamic KF
            assert(isequal(size(P0), [n n]), "ValueError: size(P0)")
            obj.P0 = P0;
            obj.Q = Q;
            assert(isequal(size(Q), [n n]))
            obj.R = R;
            assert(isequal(size(R), [ny ny]))
            obj.type = "KF";
            if nargin < 9
                obj.label = obj.type;
            end

            % Initialize variables
            obj.reset()

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Initialize estimate covariance
            obj.P = obj.P0;

            % Gain will be calculated dynamically
            obj.Kf = nan(obj.n, 1);

            % Initialize state and output estimates
            obj.xkp1_est = obj.x0;
            obj.ykp1_est = obj.C * obj.xkp1_est;

        end
        function update(obj, yk, uk)
        % obs.update(yk, uk) updates the gain and covariance matrix
        % of the Kalman filter and calculates the estimates of the
        % states and output at the current sample time.
        %
        % Arguments:
        %   yk : vector, size (ny, 1)
        %       System output measurements at current time k.
        %   uk : vector, size (nu, 1), optional
        %       System inputs at current time k.
        %

            % Update correction gain and covariance matrix
            [obj.Kf, obj.P] = kalman_update_f(obj.P, obj.A, obj.C, ...
                obj.Q, obj.R);

            % Update of state and output estimates in current timestep
            obj.xk_est = obj.xkp1_est + obj.Kf * (yk - obj.C * obj.xkp1_est);
            obj.yk_est = obj.C * obj.xk_est;

            % Prediction of state and output estimates in next timestep
            obj.xkp1_est = obj.A * obj.xk_est + obj.B * uk;
            obj.ykp1_est = obj.C * obj.xkp1_est;

        end
    end
end