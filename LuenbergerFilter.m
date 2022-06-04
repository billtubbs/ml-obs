% Kalman Filter class definition

classdef LuenbergerFilter < matlab.mixin.Copyable
% obs = LuenbergerFilter(A,B,C,D,Ts,poles,label,x0)
% Class for simulating a steady-state Kalman filter
% (i.e. with static gain).
%
% Arguments:
%   A, B, C, D : matrices
%       Discrete-time system model matrices.
%   Ts : double
%       Sampling period.
%   poles : vector
%       Desired closed-loop poles of filter dynamics.
%   label : string (optional)
%       Name.
%   x0 : vector, size(n, 1), (optional)
%       Intial state estimates.
%
% References:
%  -  D. Luenberger, "An introduction to observers," in IEEE 
%     Transactions on Automatic Control, vol. 16, no. 6, 
%     pp. 596-602, December 1971, doi: 10.1109/TAC.1971.1099826.
%
    properties (SetAccess = immutable)
        A {mustBeNumeric}
        B {mustBeNumeric}
        C {mustBeNumeric}
        D {mustBeNumeric}
        Ts {mustBeNumeric}
        poles {mustBeNumeric}
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
        function obj = LuenbergerFilter(A,B,C,D,Ts,poles,label,x0)
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.D = D;
            [n, nu, ny] = check_dimensions(A, B, C, D);
            obj.Ts = Ts;
            obj.poles = poles;
            if nargin < 7
                label = "LB";
            end
            if nargin < 8
                x0 = zeros(n, 1);
            end
            obj.label = label;
            obj.x0 = x0;
            assert(isequal(size(x0), [n 1]))
            obj.label = label;
            obj.type = "LB";

            % Compute observer gain
            if ny == 1
                obj.K = acker(A', C', poles)';
            else
                obj.K = place(A', C', poles)';
            end

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