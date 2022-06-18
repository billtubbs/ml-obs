% Abstract class definition

classdef (Abstract) AbstractLinearFilter < matlab.mixin.Copyable
% Abstract Class for inheriting when defining other
% filter classes which use a linear model, such as 
% Kalman and Luenberger filters.
%
% Properties:
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
        n (1, 1) double {mustBeInteger, mustBeNonnegative}
        nu (1, 1) double {mustBeInteger, mustBeNonnegative}
        ny (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        A double
        B double
        C double
        D double
        Ts (1, 1) double {mustBeNonnegative}
        K double
        P double
        label (1, 1) string
        x0 (:, 1) double
        xkp1_est (:, 1) double
        ykp1_est (:, 1) double
        type (1, 1) string  % is this still needed? use classdef
    end
    methods
        function obj = AbstractLinearFilter(A,B,C,D,Ts,label,x0)
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.D = D;
            [obj.n, obj.nu, obj.ny] = check_dimensions(A, B, C, D);
            obj.Ts = Ts;
            if nargin < 7
                x0 = zeros(obj.n, 1);
            else
                assert(isequal(size(x0), [obj.n 1]))
            end
            obj.type = "";
            if nargin < 6
                label = obj.type;
            end
            obj.label = label;
            obj.x0 = x0;

            % Initialize estimates
            obj.reset()

        end
        function reset(obj)
        % obj.reset()
        % Initialize variables to values specified when the filter 
        % object was created.
        %

            % Initialize state and output estimates
            obj.xkp1_est = obj.x0;
            obj.ykp1_est = obj.C * obj.xkp1_est;

        end
        function update(obj, yk, uk)
        % obs.update(yk, uk) updates the state and output estimates
        % at the next sample time.
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