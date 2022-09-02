% Steady-state Kalman Filter class definition
%
% obs = KalmanFilterSS(A,B,C,Ts,Q,R,label,x0)
% Class for simulating a steady-state Kalman filter
% (i.e. with static gain).
%
% Arguments:
%   A, B, C : matrices
%       Discrete-time system model matrices.
%   Ts : double
%       Sampling period.
%   Q : matrix, size (n, n)
%       Process noise covariance.
%   R : matrix, size (ny, ny)
%       Output measurement noise covariance.
%   label : string (optional)
%       Name.
%   x0 : vector, size(n, 1), (optional)
%       Initial state estimates.
%

classdef KalmanFilterSS < AbstractLinearFilter
    properties
        xkp1_est (:, 1) double
        ykp1_est (:, 1) double
        Pkp1 double
        K double
        Q double
        R double
    end
    methods
        function obj = KalmanFilterSS(A,B,C,Ts,Q,R,varargin)

            % Call super-class constuctor
            obj = obj@AbstractLinearFilter(A,B,C,Ts,varargin{:});
            n = obj.n;
            ny = obj.ny;

            % Set additional properties for dynamic KF
            obj.Q = Q;
            assert(isequal(size(Q), [n n]))
            obj.R = R;
            assert(isequal(size(R), [ny ny]))
            obj.type = "KFSS";
            if nargin < 7
                obj.label = obj.type;
            end

            % Compute the steady-state gain and error covariance matrix
            [obj.Pkp1,K,~,~] = idare(A',C',Q,R,[],[]);
            obj.K = K';

            % Initialize estimates
            obj.reset()

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Initialize state and output estimates
            obj.xkp1_est = obj.x0;
            obj.ykp1_est = obj.C * obj.xkp1_est;

        end
        function update(obj, yk, uk)
        % obs.update(yk, uk) updates the estimates of the
        % states and output at the next sample time.
        %
        % Arguments:
        %   yk : vector, size (ny, 1)
        %       System output measurements at current time k.
        %   uk : vector, size (nu, 1), optional
        %       System inputs at current time k.
        %

            % Update and prediction of state and output estimates
            % in next timestep
            obj.xkp1_est = obj.A * obj.xkp1_est + obj.B * uk + ...
                obj.K * (yk - obj.C * obj.xkp1_est);
            obj.ykp1_est = obj.C * obj.xkp1_est;

        end
    end
end