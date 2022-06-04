% Kalman Filter class definition

classdef KalmanFilterSS < matlab.mixin.Copyable
% obs = KalmanFilterSS(A,B,C,D,Ts,Q,R,label,x0)
% Class for simulating a steady-state Kalman filter
% (i.e. with static gain).
%
% Arguments:
%   A, B, C, D : matrices
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
%       Intial state estimates.
%
    properties (SetAccess = immutable)
        A {mustBeNumeric}
        B {mustBeNumeric}
        C {mustBeNumeric}
        D {mustBeNumeric}
        Ts {mustBeNumeric}
        Q {mustBeNumeric}
        R {mustBeNumeric}
        K {mustBeNumeric}
        P {mustBeNumeric}
        n {mustBeInteger}
        nu {mustBeInteger}
        ny {mustBeInteger}
        type
    end
    properties
        label
        x0 {mustBeNumeric}
        xkp1_est {mustBeNumeric}
        ykp1_est {mustBeNumeric}
    end
    methods
        function obj = KalmanFilterSS(A,B,C,D,Ts,Q,R,label,x0)
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.D = D;
            [n, nu, ny] = check_dimensions(A, B, C, D);
            obj.Ts = Ts;
            obj.Q = Q;
            assert(isequal(size(Q), [n n]))
            obj.R = R;
            assert(isequal(size(R), [ny ny]))
            if nargin < 8
                label = "KFSS";
            end
            if nargin < 9
                x0 = zeros(n, 1);
            end
            obj.label = label;
            obj.x0 = x0;
            assert(isequal(size(x0), [n 1]))
            obj.label = label;
            obj.type = "KFSS";

            % Model
            N = zeros(n, ny);
            G = eye(n);  % apply process noises to all states
            H = zeros(ny, n);  % no direct transmission of noises
            Gmodel = ss(A, [B G], C, [D H], Ts);

            % Use MATLAB's Kalman filter function to compute the
            % steady-state gain and covariance matrix
            [~, obj.K, obj.P] = kalman(Gmodel, Q, R, N, 'delayed');

            % Initialize estimates
            obj.xkp1_est = x0;
            obj.ykp1_est = C * obj.xkp1_est;

            % Add other useful variables
            obj.n = n;
            obj.nu = nu;
            obj.ny = ny;

        end
        function update(obj, yk, uk)
        % obs.update(yk, uk) updates the gain and covariance matrix
        % of a steady-state filter and calculates the estimates of 
        % the states and output at the next sample time.
        %
        % Arguments:
        %   yk : vector, size (ny, 1)
        %       System output measurements at current time k.
        %   uk : vector, size (nu, 1), optional
        %       System inputs at current time k.
        %

            % Update state and output estimates for next timestep
            obj.xkp1_est = obj.A * obj.xkp1_est + obj.B * uk + ...
                obj.K * (yk - obj.C * obj.xkp1_est);
            obj.ykp1_est = obj.C * obj.xkp1_est;

        end
    end
end