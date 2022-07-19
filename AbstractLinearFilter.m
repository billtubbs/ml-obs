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
%       Initial state estimates.
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
        label (1, 1) string
        x0 (:, 1) double
        type (1, 1) string  % is this still needed? use classdef
    end
    methods
        function obj = AbstractLinearFilter(A,B,C,Ts,label,x0)
            obj.A = A;
            obj.B = B;
            obj.C = C;
            [obj.n, obj.nu, obj.ny] = check_dimensions(A, B, C);
            obj.Ts = Ts;
            if nargin < 6
                x0 = zeros(obj.n, 1);
            else
                assert(isequal(size(x0), [obj.n 1]))
            end
            if nargin < 5
                label = obj.type;
            end
            obj.label = label;
            obj.x0 = x0;

        end
    end
end