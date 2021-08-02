function obs = kalman_filter(A,B,C,D,Ts,P0,Q,R,label)
% obs = kalman_filter(A,B,C,D,Ts,P0,Q,R,label)
% Creates a struct for simulating a dynamic Kalman
% filter (with time-varying gain)
%
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
    n = size(A,1);
    obs.K = nan(n,1);
    obs.static_gain = false;
    obs.xkp1_est = zeros(n,1);
    obs.ykp1_est = C * obs.xkp1_est;
end