% Simulate observers in Simulink with S-functions
%
% Simulates the same system and set of observers as in
% in test_MKF_example.m and checks the results are
% identicial

clear all

% See this Simulink model file:
sim_model = 'MKF_example_sim';

% Generate randomly-occurring shocks
% Reset random number generator
seed = 22;
rng(seed)

% Load system model
% First order SISO system with one input disturbance
sys_rodin_step

% Sequence length
nT = 100;

% Generate random shock signal
[Wp, alpha] = sample_random_shocks(nT+1, epsilon, sigma_wp(2), sigma_wp(1));

% Other inputs to system
X0 = zeros(n,1);
t = Ts*(0:nT)';
U = zeros(nT+1,1);
U(t>=5) = 1;
V = sigma_M*randn(nT+1, 1);

% Inputs to simulink model
inputs.U = [t U];
inputs.V = [t V];
inputs.Wp = [t Wp];

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Bw = B(:, ~u_meas);
Du = D(:, u_meas);

% Steady-state Kalman filters
Q = diag([0.01^2 0.1^2]);
R = 0.1^2;
obs_model = model;
obs_model.B = Bu;
obs_model.D = Du;
obs_model.Q = Q;
obs_model.R = R;
Bu = B(:,1);  % observer model without unmeasured inputs
Du = D(:,1);

% Prediction form
KFPSS = KalmanFilterPSS(obs_model,'KFPSS');

% Filtering form
KFFSS = KalmanFilterFSS(obs_model,'KFFSS');

% Dynamic Kalman filters
P0 = eye(n);

% Prediction form
KFP = KalmanFilterP(obs_model,P0,'KFP');

% Filtering form
KFF = KalmanFilterF(obs_model,P0,'KFF');

% Multi-model observer - sequence fusion
% TODO: This is not currently working
P0 = eye(n);
Q0 = diag([0.01^2 0]);
R = 0.1^2;
f = 15;  % fusion horizon
m = 1;  % maximum number of shocks
d = 3;  % spacing parameter
MKF1 = MKFObserverSF_RODD(model,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,'MKF1');

% Multi-model observer - sequence pruning
P0 = eye(n);
Q0 = diag([0.01^2 0]);
R = 0.1^2;
nh = 5;  % number of filters
n_min = 2;  % minimum life of cloned filters
MKF2 = MKFObserverSP(model,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,nh,n_min,'MKF2');

fprintf("Running Simulink simulation...\n")
sim_out = sim(sim_model, 'StopTime', string(nT*Ts), ...
    'ReturnWorkspaceOutputs', 'on');

% Check built-in Simulink Kalman filter block estimates are same as
% steady-state Matlab observer object
assert(max(abs(sim_out.X_hat_KF.Data - sim_out.X_hat_KFPSS.Data), [], [1 2]) < 1e-12)

% Check Kalman filter estimates are close to true system states
assert(mean(abs(sim_out.X_hat_KFPSS.Data - sim_out.X.Data), [1 2]) < 0.25)

% Load simulation results produced by test_MKF_example.m
% Save results - these are used by test_MKF_example_sim.m
filename = 'MKF_example.csv';
results_dir = 'results';
sim_results = readtable(fullfile(results_dir, filename));

% This loads simulation results for the following observers:
% - Xk_est_KFSS, Yk_est_KFSS
% - Xk_est_KF1, Yk_est_KF1
% - Xk_est_MKF1, Yk_est_MKF1
% - Xk_est_MKF2, Yk_est_MKF2

% Check all Simulink observer estimates are same as simulation results
% from file

% Steady-state Kalman filter - prediction form
assert(max(abs(sim_out.X_hat_KFPSS.Data - ...
    sim_results{:, {'Xk_est_KFSS_1', 'Xk_est_KFSS_2'}}), [], [1 2]) < 1e-12)

% Dynamic Kalman filter - prediction form
assert(max(abs(sim_out.X_hat_KFP.Data - sim_results{:, {'Xk_est_KF1_1', 'Xk_est_KF1_2'}}), [], [1 2]) < 1e-12)

% Multiple-model observer - SF
%assert(max(abs(sim_out.X_hat_MKF1.Data - sim_results{:, {'Xk_est_MKF1_1', 'Xk_est_MKF1_2'}}), [], [1 2]) < 1e-12)

% Multiple-model observer - SP
%TODO: Estimates seem to be the same but one step ahead
assert(max(abs(sim_out.X_hat_MKF2.Data - sim_results{:, {'Xk_est_MKF2_1', 'Xk_est_MKF2_2'}}), [], [1 2]) < 1e-12)

%disp("Simulations complete")
