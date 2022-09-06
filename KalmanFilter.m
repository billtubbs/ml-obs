% Kalman Filter class definition
%
% obs = KalmanFilter(A,B,C,D,Ts,P0,Q,R,label,x0)
% Class for simulating a dynamic Kalman filter
% (i.e. with time-varying gain and estimation error
% covariance).
%
% This is the prediction form of the KF which 
% produces prior estimates of the states and
% outputs in the next time instant given the data
% at the current time instant:
%
%  x_hat(k+1|k) : estimate of states at time k+1
%  y_hat(k+1|k) : estimate of outputs at time k+1
%
% The system model is defined as:
%
%   x(k+1) = A(k) x(k) + B(k) u(k) + w(k)
%     y(k) = C(k) x(k) + v(k)
%
% Note: there is no direct transmission (D = 0).
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

% TODO: Could rename this KalmanFilterP

classdef KalmanFilter < AbstractLinearFilter
    properties
        xkp1_est (:, 1) double
        ykp1_est (:, 1) double
        Pkp1 double
        K double
        Q double
        R double
        P0 double
    end
    methods
        function obj = KalmanFilter(A,B,C,Ts,P0,Q,R,varargin)

            % Call super-class constuctor
            obj = obj@AbstractLinearFilter(A,B,C,Ts,varargin{:});
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
            if nargin < 8
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
            obj.Pkp1 = obj.P0;

            % Gain will be calculated dynamically
            obj.K = nan(obj.n, 1);

            % Initialize state and output estimates
            obj.xkp1_est = obj.x0;
            obj.ykp1_est = obj.C * obj.xkp1_est;

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk) updates the gain and covariance matrix
        % of the Kalman filter and calculates the estimates of the
        % states and output at the next sample time given the 
        % current measurement and inputs.
        %
        % Arguments:
        %   yk : vector, size (ny, 1)
        %       System output measurements at current time k.
        %   uk : vector, size (nu, 1), optional
        %       System inputs at current time k.
        %

            % Update correction gain and covariance matrix
            [obj.K, obj.Pkp1] = kalman_update(obj.Pkp1, obj.A, obj.C, ...
                obj.Q, obj.R);

            % Update state and output estimates in next timestep
            obj.xkp1_est = obj.A * obj.xkp1_est + obj.B * uk + ...
                obj.K * (yk - obj.C * obj.xkp1_est);
            obj.ykp1_est = obj.C * obj.xkp1_est;

        end
    end
end