function LB = luenberger_filter(A,B,C,D,Ts,poles,label)
% LB = luenberger_filter(A,B,C,D,Ts,poles,label)
% Creates a struct for simulating a Luenberger
% observer (steady-state).
    LB.A = A;
    LB.B = B;
    LB.C = C;
    LB.D = D;
    LB.Ts = Ts;
    LB.poles = poles;
    LB.Gmodel = ss(A,B,C,D,Ts);
    nu = size(B, 2);
    ny = size(C, 1);
    % Compute observer gain
    if ny == 1
        LB.K = acker(A', C', poles)';
    else
        LB.K = place(A', C', poles)';
    end
    LB.label = label;
    LB.status = 1;
    n = size(A,1);
    LB.xkp1_est = zeros(n,1);
    LB.ykp1_est = C * LB.xkp1_est;
end