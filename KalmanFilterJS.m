% Kalman Filter class definition
%
% obs = KalmanFilterJS(models,P0,label,x0)
% Class for simulating a dynamic Kalman filter
% (i.e. with time-varying gain and estimation error
% covariance). This version can be used with jump
% systems where the system model switches between
% a finite set of models each time step.
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
%   models : (1, nj) cell array of structs
%       Each struct contains the parameters of a linear
%       model of the system dynamics. These include: A, B, 
%       and C for the system matrices, Q and R for the
%       state error covariance and output measurement 
%       noise covariance, and Ts for the sample period.
%   P0 : (n, n) double
%       Initial covariance matrix of the state estimates
%       (same for each filter).
%   label : string
%       Arbitrary name to identify the observer.
%   x0 : (n, 1) double
%       Initial state estimates (optional, default zeros).
%


classdef KalmanFilterJS < matlab.mixin.Copyable
    properties (SetAccess = immutable)
        nu (1, 1) double {mustBeInteger}
        ny (1, 1) double {mustBeInteger}
        n (1, 1) double {mustBeInteger}
        nj (1, 1) double {mustBeInteger}
    end
    properties
        models (1, :) cell
        K double
        P0 double
        label (1, 1) string
        x0 (:, 1) double
        rk (:, 1) double
        xkp1_est (:, 1) double
        ykp1_est (:, 1) double
        Pkp1 double
        r0 double {mustBeInteger}
        type (1, 1) string  % is this still needed? use classdef
    end
    methods
        function obj = KalmanFilterJS(models,P0,label,x0,r0)
            arguments
                models (1, :) cell
                P0 double
                label (1, 1) string = ""
                x0 = []
                r0 = 1  % default system mode at k = 0
            end

            % Get number of system models and check their dimensions
            [nj, n, nu, ny] = check_models(models);
            obj.nj = nj;
            obj.nu = nu;
            obj.n = n;
            obj.ny = ny;

            % Check dimensions of other parameters
            assert(isequal(size(P0), [n n]), "ValueError: size(P0)")
            for j = 1:nj
                assert(isequal(size(models{j}.Q), [n n]))
                assert(isequal(size(models{j}.R), [ny ny]))
            end

            % Determine initial state values
            if isempty(x0)
                x0 = zeros(obj.n, 1);  % default initial states
            else
                assert(isequal(size(x0), [obj.n 1]))
            end

            % Store parameters
            obj.models = models;
            obj.P0 = P0;
            obj.x0 = x0;
            obj.r0 = r0;
            obj.type = "KFJS";
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

            % Initialize estimate covariance
            obj.Pkp1 = obj.P0;

            % Gain will be calculated dynamically
            obj.K = nan(obj.n, 1);

            % Initialize state and output estimates
            obj.xkp1_est = obj.x0;
            obj.rk = obj.r0;
            obj.ykp1_est = obj.models{obj.rk}.C * obj.xkp1_est;

        end
        function update(obj, yk, uk, rk)
        % obj.update(yk, uk, rk) updates the gain and covariance 
        % matrix of the Kalman filter and calculates the estimates 
        % of the states and output at the next sample time given
        % the current measurement, inputs, and system mode.
        %
        % Arguments:
        %   yk : vector, size (ny, 1)
        %       System output measurements at current time k.
        %   uk : vector, size (nu, 1)
        %       System inputs at current time k.
        %   rk : vector, size (nj, 1)
        %       System mode at current time k.
        %

            % Update correction gain and covariance matrix
            [obj.K, obj.Pkp1] = kalman_update(obj.Pkp1, ...
                obj.models{rk}.A, obj.models{rk}.C, ...
                obj.models{rk}.Q, obj.models{rk}.R);

            % Update state and output estimates in next timestep
            obj.xkp1_est = obj.models{rk}.A * obj.xkp1_est ...
                + obj.models{rk}.B * uk + ...
                obj.K * (yk - obj.models{rk}.C * obj.xkp1_est);
            obj.ykp1_est = obj.models{rk}.C * obj.xkp1_est;

        end
    end
end