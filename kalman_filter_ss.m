function KF = kalman_filter_ss(A,B,C,D,Ts,Q,R,label)
% KF = kalman_filter_ss(A,B,C,D,Ts,Q,R,label)
% Creates a struct for simulating a steady-state
% Kalman filter.
    KF.A = A;
    KF.B = B;
    KF.C = C;
    KF.D = D;
    KF.Ts = Ts;
    KF.Q = Q;
    KF.R = R;
    KF.label = label;
    KF.status = 1;
    % Model
    n = size(A,1);
    ny = size(C,1);
    N = zeros(n,ny);
    Gkf = eye(n);
    Hkf = zeros(ny,n);
    Gmodel = ss(A,[B Gkf],C,[D Hkf],Ts);
    % Use MATLAB's Kalman filter object to compute
    % steady-state gain and covariance matrix
    [KF.KalmanFilter, KF.K, KF.P] = ...
        kalman(Gmodel,Q,R,N,'delayed');
    KF.xkp1_est = zeros(n,1);
    KF.ykp1_est = C * KF.xkp1_est;
end