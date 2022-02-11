% Test Simulink model and S-functions

clear all

% First run simulations in MATLAB:
test_MKF_example

% This produces simulation results in struct:
% sim_results.KF1
% sim_results.MKF1
% sim_results.MKF2

% Simulate same observers in Simulink
% See this Simulink model file:
model = 'MKF_example_sim';

% Additional parameters for simulink model
inputs.U = [t U];
inputs.V = [t V];
N = zeros(n,ny);
inputs.Wp = [t Wp];

fprintf("Running Simulink simulation...\n")
sim_out = sim(model, 'ReturnWorkspaceOutputs', 'on');

% Check both Kalman filter state estimates are the same
assert(max(abs(sim_out.X_hat_KF.Data - sim_out.X_hat_KF1.Data), [], [1 2]) < 1e-6)

% Check Kalman filter estimates are close to true system states
assert(mean(abs(sim_out.X_hat_KF.Data - sim_out.X.Data), [1 2]) < 0.5)

% Check all Simulink observer estimates are same as MATLAB estimates
assert(max(abs(sim_out.X_hat_KF1.Data - sim_results.KF1.Xk_est), [], [1 2]) < 1e-8)
assert(max(abs(sim_out.X_hat_MKF2.Data - sim_results.MKF2.Xk_est), [], [1 2]) < 1e-8)
assert(max(abs(sim_out.X_hat_MKF1.Data - sim_results.MKF1.Xk_est), [], [1 2]) < 1e-8)

% TODO: add MKF_RODD observer to simulink model

disp("Simulations complete")
