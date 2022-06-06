% Kalman Filter class definition

classdef KalmanFilter < matlab.mixin.Copyable
% obs = KalmanFilter(A,B,C,D,Ts,P0,Q,R,label,x0)
% Class for simulating a dynamic Kalman filter
% (i.e. with time-varying gain and estimation error
% covariance).
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
%       Intial state estimates.
%
    properties (SetAccess = immutable)
        Ts (1, 1) double {mustBeNonnegative}
        n (1, 1) double {mustBeInteger, mustBeNonnegative}
        nu (1, 1) double {mustBeInteger, mustBeNonnegative}
        ny (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        A double
        B double
        C double
        D double
        P0 double
        Q double
        R double
        K double
        P double
        label (1, 1) string
        x0 (:, 1) double
        xkp1_est (:, 1) double
        ykp1_est (:, 1) double
        type (1, 1) string  % is this still needed? use classdef
    end
    methods
        function obj = KalmanFilter(A,B,C,D,Ts,P0,Q,R,label,x0)
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.D = D;
            [n, nu, ny] = check_dimensions(A, B, C, D);
            obj.Ts = Ts;
            obj.P0 = P0;
            assert(isequal(size(P0), [n n]), "ValueError: size(P0)")
            obj.Q = Q;
            assert(isequal(size(Q), [n n]))
            obj.R = R;
            assert(isequal(size(R), [ny ny]))
            obj.type = "KF";
            if nargin < 10
                x0 = zeros(n, 1);
            else
                assert(isequal(size(x0), [n 1]))
            end
            obj.x0 = x0;
            if nargin < 9
                label = obj.type;
            end
            obj.label = label;

            % Add other useful variables
            obj.n = n;
            obj.nu = nu;
            obj.ny = ny;

            % Initialize estimates
            obj.reset()

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Initialize estimate covariance
            obj.P = obj.P0;

            % Initialize estimates
            obj.xkp1_est = obj.x0;
            obj.ykp1_est = obj.C * obj.xkp1_est;

            % Gain will be calculated dynamically
            obj.K = nan(obj.n, 1);

        end
        function update(obj, yk, uk)
        % obs.update(yk, uk) updates the gain and covariance matrix
        % of the Kalman filter and calculates the estimates of the
        % states and output at the next sample time.
        %
        % Arguments:
        %   yk : vector, size (ny, 1)
        %       System output measurements at current time k.
        %   uk : vector, size (nu, 1), optional
        %       System inputs at current time k.
        %

            % Update observer gain and covariance matrix
            [obj.K, obj.P] = kalman_update(obj.P, obj.A, obj.C, obj.Q, obj.R);

            % Update state and output estimates for next timestep
            obj.xkp1_est = obj.A * obj.xkp1_est + obj.B * uk + ...
                obj.K * (yk - obj.C * obj.xkp1_est);
            obj.ykp1_est = obj.C * obj.xkp1_est;

        end
    end
end