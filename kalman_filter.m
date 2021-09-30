function obs = kalman_filter(A,B,C,D,Ts,P0,Q,R,label,x0)
% obs = kalman_filter(A,B,C,D,Ts,P0,Q,R,label,x0)
% Creates a struct for simulating a dynamic Kalman
% filter (with time-varying gain)
%
% Arguments:
%	A, B, C, D : discrete-time system model matrices.
%   Ts : sample period.
%   P0 : Initial value of covariance matrix of the state
%       estimates.
%   Q : Process noise covariance matrix.
%   R : Output measurement noise covariance matrix.
%   label : string name.
%   x0 : intial state estimates (optional).
%
    n = check_dimensions(A, B, C, D);
    if nargin == 9
        % Default initial state estimate
        x0 = zeros(n, 1);
    end
    obs.A = A;
    obs.B = B;
    obs.C = C;
    obs.D = D;
    obs.Ts = Ts;
    obs.P0 = P0;
    obs.P = P0;
    obs.Q = Q;
    obs.R = R;
    obs.label = label;
    obs.status = 1;
    obs.K = nan(n, 1);
    obs.static_gain = false;
    obs.xkp1_est = x0;
    obs.ykp1_est = C * obs.xkp1_est;
    
end