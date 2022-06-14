% Kalman Filter class definition

classdef KalmanFilter < KalmanFilterSS
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
    properties
        P0 double
    end
    methods
        function obj = KalmanFilter(A,B,C,D,Ts,P0,Q,R,varargin)

            % Call super-class constuctor
            obj = obj@KalmanFilterSS(A,B,C,D,Ts,Q,R,varargin{:});

            % Set additional properties for dynamic KF
            n = obj.n;
            assert(isequal(size(P0), [n n]), "ValueError: size(P0)")
            obj.P0 = P0;
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
            obj.K = nan(obj.n, 1);

            % Initialize estimates
            reset@KalmanFilterSS(obj);

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
            update@KalmanFilterSS(obj, yk, uk);

        end
    end
end