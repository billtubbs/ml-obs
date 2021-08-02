function obs = luenberger_filter(A,B,C,D,Ts,poles,label)
% obs = luenberger_filter(A,B,C,D,Ts,poles,label)
% Creates a struct for simulating a Luenberger
% observer (steady-state).
    obs.A = A;
    obs.B = B;
    obs.C = C;
    obs.D = D;
    obs.Ts = Ts;
    obs.poles = poles;
    ny = size(C, 1);
    % Compute observer gain
    if ny == 1
        obs.K = acker(A', C', poles)';
    else
        obs.K = place(A', C', poles)';
    end
    obs.static_gain = true;
    obs.label = label;
    obs.status = 1;
    n = size(A,1);
    obs.xkp1_est = zeros(n,1);
    obs.ykp1_est = C * obs.xkp1_est;
end