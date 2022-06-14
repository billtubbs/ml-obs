% Various observers designed for the systems defined in:
%  - sys_rodin_step_2x2.m
%  - sys_rodin_step_2x2_sym.m

addpath("~/process-observers/")

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);
nw = sum(~u_meas);

% Disturbance input (used by SKF observer)
Bw = B(:, ~u_meas);
nw = sum(~u_meas);

% Check observability of system
Qobs = obsv(A,C);
unobs = length(A) - rank(Qobs);
assert(unobs == 0);

% Steady-state Luenberger observer 1
% Specify poles of observer dynamics
poles = [0.8; 0.82; 0.85; 0.85];
LB1 = luenberger_filter(A,Bu,C,Du,Ts,poles,'LB1');

% Steady-state Luenberger observer 2
% Specify poles of observer dynamics
poles = [0.6; 0.65; 0.6; 0.65];
LB2 = luenberger_filter(A,Bu,C,Du,Ts,poles,'LB2');

% Process noise covariance for states 1 and 2
% These are used by all observers
q1 = 0.01; q2 = 0.01;

% Steady-state Kalman filter 1 - tuned to input noise
Q = diag([q1 q2 sigma_wp(1,1)^2 sigma_wp(2,1)^2]);
R = diag(sigma_M.^2);
KFSS1 = KalmanFilterSS(A,Bu,C,Du,Ts,Q,R,'KFSS1');

% Steady-state Kalman filter 2 - tuned to input shocks
Q = diag([q1 q2 sigma_wp(1,2)^2 sigma_wp(2,2)^2]);
R = diag(sigma_M.^2);
KFSS2 = KalmanFilterSS(A,Bu,C,Du,Ts,Q,R,'KFSS2');

% Kalman filter 1 - tuned to sigma_wp(1)
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([q1 q2 sigma_wp(1,1)^2 sigma_wp(2,1)^2]);
R = diag(sigma_M.^2);
KF1 = KalmanFilter(A,Bu,C,Du,Ts,P0,Q,R,'KF1');

% Kalman filter 2 - tuned to sigma_wp(2)
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([q1 q2 sigma_wp(1,2)^2 sigma_wp(2,2)^2]);
R = diag(sigma_M.^2);
KF2 = KalmanFilter(A,Bu,C,Du,Ts,P0,Q,R,'KF2');

% Kalman filter 3 - manually tuned
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([q1 q2 0.1^2 0.1^2]);
R = diag(sigma_M.^2);
KF3 = KalmanFilter(A,Bu,C,Du,Ts,P0,Q,R,'KF3');

% Multiple model filter 1
label = 'MKF_SF1';
P0 = 1000*eye(n);
Q0 = diag([q1 q2 0 0]);  % TODO: Is this correct?
R = diag(sigma_M.^2);
f = 3;  % 5 fusion horizon
m = 1;  % 1 maximum number of shocks
d = 2;  % 10 spacing parameter
MKF1 = MKFObserverRODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Multiple model filter 2
label = 'MKF_SF2';
P0 = 1000*eye(n);
Q0 = diag([q1 q2 0 0]);  % TODO: Is this correct?
R = diag(sigma_M.^2);
%R = diag([1; 2.3].*sigma_M.^2);
f = 5;  % 10 fusion horizon
m = 2;  % 2 maximum number of shocks
d = 2;  % 5 spacing parameter
MKF2 = MKFObserverRODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% General MKF equivalent to MKF2
%MKF3 = mkf_filter({A,A},{B,B},{C,C},{D,D},Ts,repmat({P0},1,MKF2.n_filt),Q,R,MKF2.S,MKF2.p_seq,d,'MKF3');

% Multiple model AFMM filter 1
label = 'MKF_SP2';
P0 = 1000*eye(n);
Q0 = diag([q1 q2 0 0]);
R = diag(sigma_M.^2);
f = 100;  % sequence history length
n_filt = 10;  % number of filters
n_min = 4;  % minimum life of cloned filters
AFMM1 = MKFObserverSP(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label);

% Multiple model AFMM filter 2
label = 'MKF_SP2';
P0 = 1000*eye(n);
Q0 = diag([q1 q2 0 0]);
R = diag(sigma_M.^2);
f = 100;  % sequence history length
n_filt = 30;  % number of filters
n_min = 10;  % minimum life of cloned filters
AFMM2 = MKFObserverSP(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label);
