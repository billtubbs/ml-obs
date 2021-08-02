% Various observers designed for the system defined in:
%  - sys_rodin_step.m

% Steady-state Luenberger observer 1
% Specify poles of observer dynamics
poles = [0.8; 0.8];
LB1 = luenberger_filter(A,B,C,D,Ts,poles,'LB1');

% Steady-state Luenberger observer 2
% Specify poles of observer dynamics
poles = [0.6; 0.6];
LB2 = luenberger_filter(A,B,C,D,Ts,poles,'LB2');

% Specify covariance for state variable 1
% This is used by all observers
Q1 = 0.01;

% Adjustment factor for Q2 (applies to KF1, KF2, KF3, SKF only)
% logspace(-2, 2, 9)
% 0.0100    0.0316    0.1000    0.3162    1.0000    3.1623   10.0000   31.6228  100.0000
% logspace(-1, 1, 9)
% 0.1000    0.1778    0.3162    0.5623    1.0000    1.7783    3.1623    5.6234   10.0000
% logspace(-1, 1, 17)
% 0.5623    0.7499    1.0000    1.3335    1.7783
adj =  1;

% Steady-state Kalman filter 1 - tuned to sigma_wp(1)
Q = diag([Q1 sigma_wp(1)^2]);
R = sigma_M^2;
KFSS1 = kalman_filter_ss(A,B,C,D,Ts,Q,R,'KFSS1');

% Steady-state Kalman filter 2 - tuned to sigma_wp(2)
Q = diag([Q1 sigma_wp(2)^2]);
R = sigma_M^2;
KFSS2 = kalman_filter_ss(A,B,C,D,Ts,Q,R,'KFSS2');

% Kalman filter 1 - tuned to sigma_wp(1)
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([Q1 adj*sigma_wp(1)^2]);
R = sigma_M^2;
KF1 = kalman_filter(A,B,C,D,Ts,P0,Q,R,'KF1');

% Kalman filter 2 - tuned to sigma_wp(2)
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([Q1 adj*sigma_wp(2)^2]);
R = sigma_M^2;
KF2 = kalman_filter(A,B,C,D,Ts,P0,Q,R,'KF2');

% Kalman filter 3 - manually tuned
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([Q1 adj*0.1^2]);
R = sigma_M^2;
KF3 = kalman_filter(A,B,C,D,Ts,P0,Q,R,'KF3');

% Scheduled Kalman filter
P0 = 1000*eye(n);
Q0 = diag([Q1 1]);
R = sigma_M^2;
SKF = kalman_filter(A,B,C,D,Ts,P0,Q0,R,'SKF');
SKF.Q0 = Q0;
SKF.sigma_wp = sigma_wp;

% Multiple model filter 1
label = 'MKF1';
P0 = 1000*eye(n);
Q0 = diag([Q1 1]);
R = sigma_M^2;
f = 3;  % fusion horizon
m = 2;  % maximum number of shocks
d = 5;  % spacing parameter
MKF1 = mkf_filter_RODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Modified version of MKF1 with custom-shifted sequences
% seq_mod = cell2mat(MKF1.seq);
% seq_mod(2:4, :) = circshift(seq_mod(2:4,:),-2,2);
% seq_mod(5:6, :) = circshift(seq_mod(5:6,:),-1,2);
% seq_mod(6,2:3) = [0 1];
% MKF1m = MKF1;  % make a copy
% MKF1m.seq = mat2cell(seq_mod, ones(7,1), 9);
% MKF1m.label = strcat(MKF1.label, 'm');

% Multiple model filter 2
label = 'MKF2';
P0 = 1000*eye(n);
Q0 = diag([Q1 nan]);
R = sigma_M^2;
f = 5;  % 10 fusion horizon
m = 2;  % 2 maximum number of shocks
d = 3;  % 2 spacing parameter
MKF2 = mkf_filter_RODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Modified version of MKF1 with custom-shifted sequences
% seq_mod = cell2mat(MKF2.seq);
% seq_mod(2:11, :) = circshift(seq_mod(2:11, :), -1, 2);
% seq_mod(12:15, :) = circshift(seq_mod(12:15, :), -1, 2);
% seq_mod(16:20, 3:end) = circshift(seq_mod(16:20, 3:end), -1, 2);
% seq_mod(21:23, 3:end) = circshift(seq_mod(21:23, 3:end), -1, 2);
% seq_mod(24:28, 5:end) = circshift(seq_mod(24:28, 5:end), -1, 2);
% seq_mod(29:30, 5:end) = circshift(seq_mod(29:30, 5:end), -1, 2);
% seq_mod(31:35, 7:end) = circshift(seq_mod(31:35, 7:end), -1, 2);
% seq_mod(36:36, 7:end) = circshift(seq_mod(36:36, 7:end), -1, 2);
% seq_mod(37:41, 9:end) = circshift(seq_mod(37:41, 9:end), -1, 2);
% MKF2m = MKF2;  % make a copy
% MKF2m.seq = mat2cell(seq_mod, ones(56,1), 20);
% MKF2m.label = strcat(MKF2.label, 'm');

% General MKF equivalent to MKF2
%MKF3 = mkf_filter({A,A},{B,B},{C,C},{D,D},Ts,repmat({P0},1,MKF2.n_filt),Q,R,MKF2.S,MKF2.p_seq,d,'MKF3');
