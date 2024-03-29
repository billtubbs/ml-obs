% Steady-state Kalman Filter class definition
%
% obs = KalmanFilterPSS(model,label,x0,reset)
% Class for simulating a steady-state Kalman filter
% (i.e. with static gain). This is the prediction form of 
% the KF, which produces predictions of the states and
% outputs at the next time instant given the data at the 
% current time:
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
%   model : struct
%       Struct containing the parameters of a linear
%       model of the system dynamics. These include: A, B, 
%       and C for the system matrices, Q and R for the
%       state error covariance and output measurement 
%       noise covariance, and the sampling period, Ts.
%   label : string (optional)
%       Name.
%   x0 : vector, size(n, 1), (optional)
%       Initial state estimates at time k = 0.
%   reset : logical (default, true)
%       If true, the objects reset method is called after
%       initialization (this is mainly intended for use by
%       other objects instantiating an instance without
%       reseting).
%

classdef KalmanFilterPSS < AbstractLinearFilter
    properties
        P0 double
        K double
        Pkp1 double
    end
    methods
        function obj = KalmanFilterPSS(model,label,x0,reset)
            arguments
                model struct
                label (1, 1) string = ""
                x0 (:, 1) double = []
                reset logical = true
            end

            % Call super-class constructor
            obj = obj@AbstractLinearFilter(model,"KFPSS",label,x0,reset)

            % Check size of other parameters
            assert(isequal(size(model.Q), [obj.n obj.n]), ...
                "ValueError: size(model.Q)")
            assert(isequal(size(model.R), [obj.ny obj.ny]), ...
                "ValueError: size(model.R)")

            % Compute the steady-state gain and error covariance matrix
            % This is the gain for the filtering form of the KF:
            [Kf, obj.Pkp1] = kalman_gain_ss(obj.model.A, obj.model.C, ...
                obj.model.Q, obj.model.R);

            % The gain for the prediction form of the KF:
            obj.K = obj.model.A * Kf;

            if reset
                % Initialize variables
                obj.reset()
            end

        end
        function update(obj, yk, uk)
        % obs.update(yk, uk) updates the estimates of the
        % states and output at the next sample time:
        %         x_est(k+1|k)
        %         y_est(k+1|k)
        % given the current measurement y(k) and the current 
        % input u(k).
        %
        % Arguments:
        %   yk : vector, size (ny, 1)
        %       System output measurements at current time k.
        %   uk : vector, size (nu, 1), optional
        %       System inputs at current time k.
        %
        % After calling this method, the following properties will
        % be updated:
        %   - obs.xkp1_est  (x_est(k+1|k))
        %   - obs.ykp1_est  (y_est(k+1|k))
        %

            % Update predictions of states and outputs in 
            % next timestep
            obj.xkp1_est = obj.model.A * obj.xkp1_est ...
                + obj.model.B * uk ...
                + obj.K * (yk - obj.model.C * obj.xkp1_est);
            obj.ykp1_est = obj.model.C * obj.xkp1_est;

        end
    end
end