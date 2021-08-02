function KF = kalman_filter(A,B,C,D,Ts,P0,Q,R,label)
% KF = kalman_filter(A,B,C,D,Ts,P0,Q,R,label)
% Creates a struct for simulating a time-varying
% Kalman filter.
    KF.A = A;
    KF.B = B;
    KF.C = C;
    KF.D = D;
    KF.Ts = Ts;
    KF.P0 = P0;
    KF.P = P0;
    KF.Q = Q;
    KF.R = R;
    KF.label = label;
    KF.status = 1;
    n = size(A,1);
    KF.K = nan(n,1);
    KF.xkp1_est = zeros(n,1);
    KF.ykp1_est = C * KF.xkp1_est;
end