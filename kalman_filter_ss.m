function obs = kalman_filter_ss(A,B,C,D,Ts,Q,R,label,x0)
% obs = kalman_filter_ss(A,B,C,D,Ts,Q,R,label,x0)
% Creates a struct for simulating a steady-state
% Kalman filter (i.e. with static gain).
%
% Arguments:
%   A, B, C, D : discrete-time system model matrices.
%   Ts : sample period.
%   Q : Process noise covariance matrix.
%   R : Output measurement noise covariance matrix.
%   label : string name.
%   x0 : intial state estimates (optional).
%
    [n, nu, ny] = check_dimensions(A, B, C, D);
    if nargin == 8
        x0 = zeros(n, 1);
    end
    assert(isequal(size(x0), [n 1]))
    obs.A = A;
    obs.B = B;
    obs.C = C;
    obs.D = D;
    obs.Ts = Ts;
    obs.Q = Q;
    assert(isequal(size(Q), [n n]))
    obs.R = R;
    assert(isequal(size(R), [ny ny]))
    obs.label = label;
    obs.status = 1;
    obs.type = "KFSS";

    % Model
    N = zeros(n, ny);
    G = eye(n);  % apply process noises to all states
    H = zeros(ny, n);  % no direct transmission of noises
    Gmodel = ss(A, [B G], C, [D H], Ts);

    % Use MATLAB's Kalman filter object to compute
    % steady-state gain and covariance matrix
    [obs.sys, obs.K, obs.P] = ...
        kalman(Gmodel, Q, R, N, 'delayed');

    % Initialize estimates
    obs.xkp1_est = x0;
    obs.ykp1_est = C * obs.xkp1_est;

    % Flag used buy update_KF function
    obs.static_gain = true;

    % Add other useful variables
    obs.n = n;
    obs.nu = nu;
    obs.ny = ny;

end