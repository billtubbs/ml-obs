% Various observers designed for the system defined in:
%  - sys_rodin_step.m

addpath("~/process-observers/")

% Check observability of system
Qobs = obsv(A, C);
unobs = length(A) - rank(Qobs);
assert(unobs == 0);

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Bw = B(:, ~u_meas);
Du = D(:, u_meas);

% Steady-state Luenberger observer 1
% Specify poles of observer dynamics
poles = [0.8; 0.8];
LB1 = LuenbergerFilter(A,Bu,C,Du,Ts,poles,'LB1');

% Steady-state Luenberger observer 2
% Specify poles of observer dynamics
poles = [0.6; 0.6];
LB2 = LuenbergerFilter(A,Bu,C,Du,Ts,poles,'LB2');

% Specify covariance for state variable 1
% This is used by all observers
q1 = 0.01;

% Steady-state Kalman filter 1 - tuned to sigma_wp(1)
Q = diag([q1 sigma_wp(1)^2]);
R = sigma_M^2;
KFSS1 = KalmanFilterSS(A,Bu,C,Du,Ts,Q,R,'KFSS1');

% Steady-state Kalman filter 2 - tuned to sigma_wp(2)
Q = diag([q1 sigma_wp(2)^2]);
R = sigma_M^2;
KFSS2 = KalmanFilterSS(A,Bu,C,Du,Ts,Q,R,'KFSS2');

% Kalman filter 1 - tuned to sigma_wp(1)
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([q1 sigma_wp(1)^2]);
R = sigma_M^2;
KF1 = KalmanFilter(A,Bu,C,Du,Ts,P0,Q,R,'KF1');

% Kalman filter 2 - tuned to sigma_wp(2)
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([q1 sigma_wp(2)^2]);
R = sigma_M^2;
KF2 = KalmanFilter(A,Bu,C,Du,Ts,P0,Q,R,'KF2');

% Kalman filter 3 - manually tuned
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([q1 0.1^2]);
R = sigma_M^2;
KF3 = KalmanFilter(A,Bu,C,Du,Ts,P0,Q,R,'KF3');

% Multiple model observer with sequence fusion #1
label = 'MKF_SF1';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
f = 3;  % fusion horizon
m = 2;  % maximum number of shocks
d = 5;  % spacing parameter
MKF_SF1 = MKFObserverSF(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Multiple model observer with sequence fusion #2
label = 'MKF_SF2';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
f = 10;  % fusion horizon
m = 1;  % maximum number of shocks
d = 1;  % spacing parameter
MKF_SF2 = MKFObserverSF(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Multiple model observer with sequence fusion based on method
% described in Robertson et al. 1995.
label = 'MKF_SF95';
MKF_SF95 = MKFObserverSF95(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% General MKF - should be equivalent to MKF_SF2
Q1 = diag([q1 sigma_wp(1)^2]);
Q2 = diag([q1 sigma_wp(2)^2]);
MKF3 = MKFObserver({A,A},{B,B},{C,C},{D,D},Ts,P0,{Q1,Q2},{R,R}, ...
    MKF_SF2.seq,MKF_SF2.T,'MKF3');
% TODO: Allow P0 to be replaced with repmat({P0},1,MKF2.n_filt)

% Multiple model observer with sequence pruning #1
label = 'MKF_SP1';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
f = 100;  % sequence history length
n_filt = 7;  % number of filters
n_min = 3;  % minimum life of cloned filters
MKF_SP1 = MKFObserverSP(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label);

% Multiple model observer with sequence pruning #2
label = 'MKF_SP2';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
f = 100;  % sequence history length
n_filt = 10;  % number of filters
n_min = 5;  % minimum life of cloned filters
MKF_SP2 = MKFObserverSP(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label);

observers = {LB1, LB2, KFSS1, KFSS2, KF1, KF2, KF3, MKF_SF1, MKF_SF2, ...
    MKF3, MKF_SP1, MKF_SP2};
