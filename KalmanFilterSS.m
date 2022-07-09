% Steady-state Kalman Filter class definition

classdef KalmanFilterSS < AbstractLinearFilter
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
    properties
        Q double
        R double
    end
    methods
        function obj = KalmanFilterSS(A,B,C,D,Ts,Q,R,varargin)

            % Call super-class constuctor
            obj = obj@AbstractLinearFilter(A,B,C,D,Ts,varargin{:});
            n = obj.n;
            ny = obj.ny;

            % Set additional properties for dynamic KF
            obj.Q = Q;
            assert(isequal(size(Q), [n n]))
            obj.R = R;
            assert(isequal(size(R), [ny ny]))
            obj.type = "KFSS";
            if nargin < 8
                obj.label = obj.type;
            end

            % Compute the steady-state gain and error covariance matrix
            [obj.P,K,~,~] = idare(A',C',Q,R,[],[]);
            obj.K = K';

            % Initialize estimates
            obj.reset()

        end
    end
end