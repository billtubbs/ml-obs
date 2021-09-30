function obs = kalman_filter_ss(A,B,C,D,Ts,Q,R,label,x0)
% obs = kalman_filter_ss(A,B,C,D,Ts,Q,R,label,x0)
% Creates a struct for simulating a steady-state
% Kalman filter (i.e. with static gain).
%
% Arguments:
%	A, B, C, D : discrete-time system model matrices.
%   Ts : sample period.
%   Q : Process noise covariance matrix.
%   R : Output measurement noise covariance matrix.
%   label : string name.
%   x0 : intial state estimates (optional).
%
    n = size(A,1);
    ny = size(C,1);
    if nargin == 8
        x0 = zeros(n,1);
    end
    obs.A = A;
    obs.B = B;
    obs.C = C;
    obs.D = D;
    obs.Ts = Ts;
    obs.Q = Q;
    obs.R = R;
    obs.label = label;
    obs.status = 1;
    % Model
    N = zeros(n,ny);
    Gkf = eye(n);
    Hkf = zeros(ny,n);
    Gmodel = ss(A,[B Gkf],C,[D Hkf],Ts);
    % Use MATLAB's Kalman filter object to compute
    % steady-state gain and covariance matrix
    [obs.KalmanFilter, obs.K, obs.P] = ...
        kalman(Gmodel,Q,R,N,'delayed');
    obs.static_gain = true;
    obs.xkp1_est = x0;
    obs.ykp1_est = C * obs.xkp1_est;
end