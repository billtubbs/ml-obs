% Various observers designed for the system defined in:
%  - sys_rodin_step_2x2.m

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Du = D(:, u_meas);

% Check observability of system
Qobs = obsv(A,C);
unobs = length(A) - rank(Qobs);
assert(unobs == 0);

% Steady-state Luenberger observer 1
% Specify poles of observer dynamics
poles = [0.8; 0.82; 0.85; 0.85];
LB1 = luenberger_filter(A,B,C,D,Ts,poles,'LB1');

% Steady-state Luenberger observer 2
% Specify poles of observer dynamics
poles = [0.6; 0.65; 0.6; 0.65];
LB2 = luenberger_filter(A,B,C,D,Ts,poles,'LB2');

% Adjustment factor for Q2 (applies to KF1, KF2, KF3, SKF, MKF1, MKF2)
% logspace(-2, 2, 9)
% 0.0100    0.0316    0.1000    0.3162    1.0000    3.1623   10.0000   31.6228  100.0000
% logspace(-1, 1, 9)
% 0.1000    0.1778    0.3162    0.5623    1.0000    1.7783    3.1623    5.6234   10.0000
% logspace(-1, 1, 17)
% 0.4217    0.5623    0.7499    1.0000    1.3335    1.7783    2.3714    3.1623    4.2170    5.6234

% For adjusting Q(3,3) and Q(4,4) on the MKF and AFMM filters
adj = 1;

% For adjusting the measurement noise covariances on all filters
Radj = 1;

% Process noise covariance for states 1 and 2
% These are used by all observers
Q1 = 0.01; Q2 = 0.01;

% Steady-state Kalman filter 1 - tuned to input noise
Q = diag([Q1 Q2 sigma_wp(1,1)^2 sigma_wp(2,1)^2]);
R = diag(sigma_M.^2);
KFSS1 = kalman_filter_ss(A,B,C,D,Ts,Q,R,'KFSS1');

% Steady-state Kalman filter 2 - tuned to input shocks
Q = diag([Q1 Q2 sigma_wp(1,2)^2 sigma_wp(2,2)^2]);
R = diag(sigma_M.^2);
KFSS2 = kalman_filter_ss(A,B,C,D,Ts,Q,R,'KFSS2');

% Kalman filter 1 - tuned to sigma_wp(1)
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([Q1 Q2 adj*sigma_wp(1,1)^2 adj*sigma_wp(2,1)^2]);
R = diag(sigma_M.^2);
KF1 = kalman_filter(A,B,C,D,Ts,P0,Q,R,'KF1');

% Kalman filter 2 - tuned to sigma_wp(2)
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([Q1 Q2 adj*sigma_wp(1,2)^2 adj*sigma_wp(2,2)^2]);
R = diag(sigma_M.^2);
KF2 = kalman_filter(A,B,C,D,Ts,P0,Q,R,'KF2');

% Kalman filter 3 - manually tuned
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([Q1 Q2 0.1^2 0.1^2]);
R = diag(Radj*sigma_M.^2);
KF3 = kalman_filter(A,B,C,D,Ts,P0,Q,R,'KF3');

% Scheduled Kalman filter
P0 = 1000*eye(n);
Q0 = diag([Q1 Q2 nan nan]);
R = diag(Radj*sigma_M.^2);
SKF = kalman_filter(A,B,C,D,Ts,P0,Q0,R,'SKF');
SKF.Q0 = Q0;
SKF.sigma_wp = sigma_wp;

% Multiple model filter 1
label = 'MKF1';
P0 = 1000*eye(n);
Q0 = diag([Q1 Q2 adj adj]);  % TODO: Is this correct?
R = diag(Radj*sigma_M.^2);
f = 3;  % 5 fusion horizon
m = 1;  % 1 maximum number of shocks
d = 2;  % 10 spacing parameter
MKF1 = mkf_filter_RODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Multiple model filter 2
label = 'MKF2';
P0 = 1000*eye(n);
Q0 = diag([Q1 Q2 adj adj]);  % TODO: Is this correct?
R = diag(Radj*sigma_M.^2);
%R = diag([1; 2.3].*sigma_M.^2);
f = 5;  % 10 fusion horizon
m = 2;  % 2 maximum number of shocks
d = 2;  % 5 spacing parameter
MKF2 = mkf_filter_RODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% General MKF equivalent to MKF2
%MKF3 = mkf_filter({A,A},{B,B},{C,C},{D,D},Ts,repmat({P0},1,MKF2.n_filt),Q,R,MKF2.S,MKF2.p_seq,d,'MKF3');

% Multiple model AFMM filter 1
label = 'AFMM1';
P0 = 1000*eye(n);
Q0 = diag([Q1 Q2 adj adj]);
R = diag(Radj*sigma_M.^2);
f = 100;  % sequence history length
n_filt = 10;  % number of filters
n_min = 3;  % minimum life of cloned filters
AFMM1 = mkf_filter_AFMM(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label);

% Multiple model AFMM filter 2
label = 'AFMM2';
P0 = 1000*eye(n);
Q0 = diag([Q1 Q2 adj adj]);
R = diag(Radj*sigma_M.^2);
f = 100;  % sequence history length
n_filt = 30;  % number of filters
n_min = 10;  % minimum life of cloned filters
AFMM2 = mkf_filter_AFMM(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label);