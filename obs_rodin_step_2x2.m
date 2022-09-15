% Various observers designed for the systems defined in:
%  - sys_rodin_step_2x2.m
%  - sys_rodin_step_2x2_sym.m

addpath("~/process-observers/")

assert(exist("A", 'var'), strcat("System model not defined. ", ...
    "Run script 'sys_rodin_step_2x2sym.m' first."))

% Check observability of system
Qobs = obsv(A,C);
unobs = length(A) - rank(Qobs);
assert(unobs == 0);

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);
nw = sum(~u_meas);

% Disturbance input (used by SKF observer)
Bw = B(:, ~u_meas);
nw = sum(~u_meas);

% Steady-state Luenberger observer 1
% Specify poles of observer dynamics
poles = [0.8; 0.82; 0.85; 0.85];
LB1 = LuenbergerFilter(A,Bu,C,Ts,poles,'LB1');

% Steady-state Luenberger observer 2
% Specify poles of observer dynamics
poles = [0.6; 0.65; 0.6; 0.65];
LB2 = LuenbergerFilter(A,Bu,C,Ts,poles,'LB2');

% Process noise covariance for states 1 and 2
% These are used by all observers
q1 = 0.01; q2 = 0.01;

% Different values for covariance matrix
Q1 = diag([q1 q2 sigma_wp(1,1)^2 sigma_wp(2,1)^2]);
Q2 = diag([q1 q2 sigma_wp(1,2)^2 sigma_wp(2,2)^2]);
Q3 = diag([q1 q2 0.0075 0.0075]);

% Covariance of output errors
R = diag(sigma_M.^2);

% Steady-state Kalman filter 1 - tuned to input noise
KFSS1 = KalmanFilterSS(A,Bu,C,Ts,Q1,R,'KFSS1');

% Steady-state Kalman filter 2 - tuned to input shocks
KFSS2 = KalmanFilterSS(A,Bu,C,Ts,Q2,R,'KFSS2');

% Observer models for new observer functions
models = {struct, struct};
models{1}.A = A;
models{1}.B = Bu;
models{1}.C = C;
models{1}.Ts = Ts;
models{1}.Q = Q1;
models{1}.R = R;
models{2}.A = A;
models{2}.B = Bu;
models{2}.C = C;
models{2}.Ts = Ts;
models{2}.Q = Q2;
models{2}.R = R;

% Parameters for manually-tuned Kalman filter (KF3)
model3 = models{1};  % makes copy
model3.Q = Q3;

% Kalman filter 1 - tuned to sigma_wp(1)
% Covariance matrices
P0 = 1000*eye(n);
KF1 = KalmanFilterF(models{1},P0,'KF1');

% Kalman filter 2 - tuned to sigma_wp(2)
% Covariance matrices
P0 = 1000*eye(n);
KF2 = KalmanFilterF(models{2},P0,'KF2');

% Kalman filter 3 - manually tuned
% Covariance matrices
P0 = 1000*eye(n);
KF3 = KalmanFilterF(model3,P0,'KF3');

% % Multiple model observer with sequence fusion #1
% label = 'MKF_SF1';
% P0 = 1000*eye(n);
% Q0 = diag([q1 q2 0 0]);  % TODO: Is this correct?
% R = diag(sigma_M.^2);
% f = 25;  % fusion horizon
% m = 2;  % maximum number of shocks
% d = 5;  % spacing parameter
% MKF_SF1 = MKFObserverSF(A,B,C,Ts,u_meas,P0,epsilon,sigma_wp, ...
%     Q0,R,f,m,d,label);
% 
% % Multiple model observer with sequence fusion #2
% label = 'MKF_SF2';
% P0 = 1000*eye(n);
% Q0 = diag([q1 q2 0 0]);  % TODO: Is this correct?
% R = diag(sigma_M.^2);
% %R = diag([1; 2.3].*sigma_M.^2);
% f = 30;  % fusion horizon
% m = 2;  % maximum number of shocks
% d = 10;  % spacing parameter
% MKF_SF2 = MKFObserverSF(A,B,C,Ts,u_meas,P0,epsilon,sigma_wp, ...
%     Q0,R,f,m,d,label);
% 
% % Multiple model observer with sequence fusion based on method
% % described in Robertson et al. 1995.
% label = 'MKF_SF95';
% MKF_SF95 = MKFObserverSF95(A,B,C,Ts,u_meas,P0,epsilon,sigma_wp, ...
%     Q0,R,f,m,d,label);
% 
% General MKF - should be equivalent to MKF2
% Q1 = diag([q1 q2 sigma_wp(1,1)^2 sigma_wp(2,1)^2]);
% Qnj = diag([q1 q2 sigma_wp(1,2)^2 sigma_wp(2,2)^2]);
% MKF3 = MKFObserverS({A,A},{B,B},{C,C},Ts,P0,{Q1,Qnj},{R,R}, ...
%     MKF_SF2.seq,MKF_SF2.T,'MKF3');
% TODO: Allow P0 to be replaced with repmat({P0},1,MKF2.n_filt)

% Multiple model observer with sequence pruning #1
label = 'MKF_SP1';
P0 = 1000*eye(n);
Q0 = diag([q1 q2 0 0]);
R = diag(sigma_M.^2);
nh = 19;  % number of filters
n_min = 5;  % minimum life of cloned filters
MKF_SP1 = MKFObserverSP(model,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,nh,n_min,label);

% Multiple model observer with sequence pruning #2
label = 'MKF_SP2';
P0 = 1000*eye(n);
Q0 = diag([q1 q2 0 0]);
R = diag(sigma_M.^2);
nh = 25;  % number of filters
n_min = 9;  % minimum life of cloned filters
MKF_SP2 = MKFObserverSP(model,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,nh,n_min,label);

% TODO: Restore
% observers = {LB1, LB2, KFSS1, KFSS2, KF1, KF2, KF3, MKF_SF1, MKF_SF2, ...
%    MKF3, MKF_SP1, MKF_SP2};
observers = {LB1, LB2, KFSS1, KFSS2, KF1, KF2, KF3, MKF_SP1, MKF_SP2};
