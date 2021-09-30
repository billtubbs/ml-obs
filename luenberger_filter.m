function obs = luenberger_filter(A,B,C,D,Ts,poles,label,x0)
% obs = luenberger_filter(A,B,C,D,Ts,poles,label,x0)
% Creates a struct for simulating a Luenberger
% observer (steady-state).
%
% Arguments:
%	A, B, C, D : discrete-time system model matrices.
%   Ts : sample period.
%   poles : Desired closed-loop poles of filter dynamics.
%   label : string name.
%   x0 : intial state estimates (optional).
%
    n = size(A,1);
    if nargin == 7
        x0 = zeros(n,1);
    end
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
    obs.xkp1_est = x0;
    obs.ykp1_est = C * obs.xkp1_est;
end