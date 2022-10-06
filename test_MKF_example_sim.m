% Simulate observers in Simulink with S-functions
%
% Simulates the same system and set of observers as in
% in test_MKF_example.m and checks the results are
% identicial

clear all

%show_plots = true;
show_plots = false;

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

% Simulate system in MATLAB
[Y,T,X] = lsim(Gpss,[U Wp],t,X0);
Ym = Y + V;  % measurement

if show_plots
    figure(2)
    subplot(2,1,1)
    plot(t,Y,t,Ym); grid on
    ylabel('y(k) and y_m(k)')
    legend('y(k)', 'ym(k)')
    title('Outputs')
    subplot(2,1,2)
    stairs(t, [U Wp]); grid on
    xlabel('k')
    ylabel('u(k) and wp(k)')
    legend('u(k)', 'wp(k)')
    title('Inputs')
end

% Simulate system in MATLAB
[Y,T,X] = lsim(Gpss,[U Wp],t,X0);
Ym = Y + V;  % measurement

if show_plots
    figure(2)
    subplot(2,1,1)
    plot(t,Y,t,Ym); grid on
    ylabel('y(k) and y_m(k)')
    legend('y(k)', 'ym(k)')
    title('Outputs')
    subplot(2,1,2)
    stairs(t, [U Wp]); grid on
    xlabel('k')
    ylabel('u(k) and wp(k)')
    legend('u(k)', 'wp(k)')
    title('Inputs')
end

% Inputs to simulink model
inputs.U = [t U];
inputs.V = [t V];
inputs.Wp = [t Wp];

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Bw = B(:, ~u_meas);
Du = D(:, u_meas);

% Kalman filter parameters
Q = diag([0.01^2 0.1^2]);
R = 0.1^2;
obs_model = model;
obs_model.B = Bu;
obs_model.D = Du;
obs_model.Q = Q;
obs_model.R = R;

% Prediction form
KFPSS = KalmanFilterPSS(obs_model,'KFPSS');

% Filtering form
KFFSS = KalmanFilterFSS(obs_model,'KFFSS');

% Dynamic Kalman filters
P0 = eye(n);

% Prediction form
KFP = KalmanFilterP(obs_model,P0,'KFP');
KF_old = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,'KF_old');

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
    sim_results{:, {'Xkp1_est_KFPSS_1', 'Xkp1_est_KFPSS_2'}}), [], [1 2]) < 1e-12)

% TODO: Remove this
observers = {KF_old};
[Y,t,X] = lsim(Gpss,[U Wp],t);
assert(max(abs(X - sim_out.X.Data), [], [1 2]) < 1e-12)
Ym = Y + V;
assert(max(abs(Ym - sim_out.Y.Data)) < 1e-12)
[Xk_est_old,Yk_est_old,DiagPk,MKF_vars] = ...
      run_simulation_obs(Ym,U,alpha,[],observers,[]);

% % Check estimates of old KF struct same as data on file
% assert(max(abs(Xk_est_old - ...
%     sim_results{:, {'Xkp1_est_KFP_1', 'Xkp1_est_KFP_2'}}), [], [1 2]) < 1e-12)

% Check KFP is the same as old KF struct
assert(max(abs(sim_out.X_hat_KFP.Data - Xk_est_old), [], [1 2]) < 1e-12)

% Dynamic Kalman filter - prediction form
assert(max(abs(sim_out.X_hat_KFP.Data - ...
    sim_results{:, {'Xkp1_est_KFP_1', 'Xkp1_est_KFP_2'}}), [], [1 2]) < 1e-12)

% Multiple-model observer - SF
%assert(max(abs(sim_out.X_hat_MKF1.Data - sim_results{:, {'Xk_est_MKF1_1', 'Xk_est_MKF1_2'}}), [], [1 2]) < 1e-12)

% Multiple-model observer - SP
%TODO: Estimates seem to be the same but one step ahead
assert(max(abs(sim_out.X_hat_MKF2.Data - sim_results{:, {'Xk_est_MKF2_1', 'Xk_est_MKF2_2'}}), [], [1 2]) < 1e-12)

%disp("Simulations complete")
