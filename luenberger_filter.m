function obs = luenberger_filter(A,B,C,D,Ts,poles,label,x0)
% obs = luenberger_filter(A,B,C,D,Ts,poles,label,x0)
% Creates a struct for simulating a Luenberger
% observer (steady-state).
%
% Arguments:
%   A, B, C, D : discrete-time system model matrices.
%   Ts : sample period.
%   poles : Desired closed-loop poles of filter dynamics.
%   label : string name.
%   x0 : intial state estimates (optional).
%
% References:
%  -  D. Luenberger, "An introduction to observers," in IEEE 
%     Transactions on Automatic Control, vol. 16, no. 6, 
%     pp. 596-602, December 1971, doi: 10.1109/TAC.1971.1099826.
%
    [n, nu, ny] = check_dimensions(A, B, C, D);
    if nargin == 7
        x0 = zeros(n,1);
    end
    obs.A = A;
    obs.B = B;
    obs.C = C;
    obs.D = D;
    obs.Ts = Ts;
    obs.poles = poles;
    obs.label = label;
    obs.status = 1;
    obs.type = "LB";

    % Compute observer gain
    if ny == 1
        obs.K = acker(A', C', poles)';
    else
        obs.K = place(A', C', poles)';
    end

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