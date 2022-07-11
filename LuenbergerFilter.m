% Luenberger filter class definition
%
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
%       Initial state estimates.
%
% References:
%  -  D. Luenberger, "An introduction to observers," in IEEE 
%     Transactions on Automatic Control, vol. 16, no. 6, 
%     pp. 596-602, December 1971, doi: 10.1109/TAC.1971.1099826.
%

classdef LuenbergerFilter < AbstractLinearFilter
    properties
        xkp1_est (:, 1) double
        ykp1_est (:, 1) double
        K double
        poles {mustBeNumeric}
    end
    methods
        function obj = LuenbergerFilter(A,B,C,D,Ts,poles,varargin)

            % Call super-class constuctor
            obj = obj@AbstractLinearFilter(A,B,C,D,Ts,varargin{:});
            n = obj.n;
            ny = obj.ny;

            % Set properties for Luenberger filter
            obj.poles = poles;
            obj.type = "LB";
            if nargin < 7
                obj.label = obj.type;
            end

            % Compute observer gain
            if ny == 1
                obj.K = acker(A', C', poles)';
            else
                obj.K = place(A', C', poles)';
            end

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
