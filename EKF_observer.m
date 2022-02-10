function obs = EKF_observer(n,f,h,params,u_meas,y_meas,dfdx,dhdx,Ts, ...
    P0,Q,R,label,x0,y0)
% obs = EKF_observer(n,f,h,params,u_meas,y_meas,dfdx,dhdx,Ts, ...
%     P0,Q,R,label,x0,y0)
% Creates a struct for simulating a discrete-time
% extended Kalman filter (EKF) for online state estimation
% of non-linear systems.
%
% System representation:
%
% x(k+1) = f(x(k), u(k), ...) + w(k)
%   y(k) = h(x(k), ...) + v(k)
%
% where w(k) and v(k) are process noise and measurement
% noise, respectively, and are zero-mean, normally-
% distributed white noises with covariances Q and R.
%
% Arguments:
%   n : Number of model states.
%   f : State transition function (non-linear).
%   h : Measurement function.
%   params : cell array containing any additional 
%       parameters that should be passed to functions f, 
%       h, dfdx, and dhdx as the final arguments.
%   u_meas : array indicating which inputs are measured.
%   y_meas : array indicating which outputs are measured.
%   dfdx : function to generate the Jacobian matrix of 
%       the state transition function.
%   dhdx : function to generate the Jacobian matrix of
%       the measurement function.
%   Ts : sample period.
%   P0 : Initial value of covariance matrix of the state
%       estimates.
%   Q : Process noise covariance matrix.
%   R : Output measurement noise covariance matrix.
%   label : string name.
%   x0 : intial state estimates (optional).
%   y0 : intial output estimates (optional).
%
    obs.n = n;
    if nargin == 13
        % Default initial state estimate
        x0 = zeros(n, 1);
    end
    if nargin < 15
        % Default initial state estimate
        y0 = zeros(size(R, 1), 1);
    end
    assert(isequal(size(x0), [n 1]), "ValueError: size(x0)")
    ny = size(y0, 1);
    obs.f = f;
    obs.h = h;
    obs.params = params;
    assert(isa(obs.f, 'function_handle'), "ValueError: f")
    assert(isa(obs.h, 'function_handle'), "ValueError: h")
    obs.u_meas = u_meas;  % TODO implement these
    obs.y_meas = y_meas;  % 
    obs.dfdx = dfdx;
    obs.dhdx = dhdx;
    assert(isa(obs.dfdx, 'function_handle'), "ValueError: dfdx")
    assert(isa(obs.dhdx, 'function_handle'), "ValueError: dhdx")
    obs.Ts = Ts;
    obs.P0 = P0;
    assert(isequal(size(P0), [n n]), "ValueError: size(P0)")
    obs.Q = Q;
    assert(isequal(size(Q), [n n]), "ValueError: size(Q)")
    obs.R = R;
    assert(isequal(size(R), [ny ny]), "ValueError: size(R)")
    obs.label = label;
    obs.status = 1;
    obs.K = nan(n,1);
    obs.xkp1_est = x0;
    obs.ykp1_est = y0;
    obs.Pkp1 = P0;
    obs.static_gain = false;
    obs.xk_est = nan(size(x0));
    obs.yk_est = nan(size(y0));
    obs.P = nan(size(P0));
    obs.type = "EKF";

end