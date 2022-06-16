% Kalman Filter class definition

classdef LuenbergerFilter < AbstractLinearFilter
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
    properties
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
    end
end
