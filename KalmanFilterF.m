% Kalman Filter class definition
%
% obs = KalmanFilterF(A,B,C,D,Ts,P0,Q,R,label,x0)
% Class for simulating a dynamic Kalman filter
% (i.e. with time-varying gain and estimation error
% covariance). This is the filtering form of the
% KF, which produces posterior estimates of the states and
% outputs at the current time instant given the data
% at the current time:
%
%  x_hat(k|k) : estimate of states at time k
%  y_hat(k|k) : estimate of outputs at time k
%
% Prior estimates of the states and outputs at the next
% time instant given the data at the current time are
% also calculated:
%
%  xkp1_hat(k+1|k) : estimate of states at time k + 1
%  ykp1_hat(k+1|k) : estimate of outputs at time k + 1
%
% Arguments:
%   A, B, C, D : matrices
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

classdef KalmanFilterF < AbstractLinearFilter
    properties
        xk_est (:, 1) double
        yk_est (:, 1) double
        xkp1_est (:, 1) double
        ykp1_est (:, 1) double
        Kf double
        Pk double
        Pkp1 double
        Sk double
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
            obj.type = "KFF";
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

            % Initialize state and output estimates
            % Note: At initialization time k = 0, xkp1_est and
            % ykp1_est represent prior estimates of the states,
            % i.e. x_est(k|k-1) and y_est(k|k-1).
            obj.xkp1_est = obj.x0; 
            obj.ykp1_est = obj.C * obj.xkp1_est;

            % Initialize error covariance P(k|k-1)
            obj.Pkp1 = obj.P0;

            % At initialization time k = 0, x_est(k|k)
            % and y_est(k|k) have not yet been computed.
            obj.xk_est = nan(obj.n, 1);
            obj.yk_est = nan(obj.ny, 1);
            obj.Pk = nan(obj.n);
            obj.Sk = nan(obj.ny);

            % Gain will be calculated at first update
            obj.Kf = nan(obj.n, 1);

        end
        function update(obj, yk, uk)
        % obs.update(yk, uk) updates the gain and covariance matrix
        % of the Kalman filter and calculates the estimates of the
        % states and output at the current sample time.
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
            [obj.Kf, obj.xk_est, obj.Pk, obj.yk_est, obj.Sk] = ...
                kalman_update_f(obj.C, obj.R, obj.xkp1_est, obj.Pkp1, yk);

            % Make predictions at next time instant
            [obj.xkp1_est,obj.ykp1_est,obj.Pkp1] = ...
                kalman_predict_f(obj.A,obj.B,obj.C,obj.Q,obj.xk_est, ...
                    obj.Pk,uk);

        end
    end
end