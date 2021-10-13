% Test EKF_observer.m and update_EKF.m
%
% To run this test use the following command:
%
% >> runtests test_EKF_observer
%

clear all


%% Example from homework from GEL-7029 course

addpath('~/process-models/pend/')

filename = 'pend_sim_benchmark.csv';
results_dir = 'results';
bench_sim_results = readtable(fullfile(results_dir, filename), ...
    'PreserveVariableNames',true);
% In MATLAB 2020b 'PreserveVariableNames' is deprecated.  Use
% readtable(filename,'PreserveVariableNames',true)

% Simulation parameters
N = 150;  % no. of simulation points = 15 sec
Ts = 0.1;  % sampling period

% Pendulum parameters
params.K = 1.2;
params.m = 0.3;
params.L = 0.4;
params.g = 9.8;
params.dt = Ts;

% Eqns for augmented model with one integrator
na = 3;
f = @pend_StateFcn;
h = @pend_MeasurementFcn;
u_meas = true;
y_meas = true;
dfdx = @pend_F;
dhdx = @pend_H;
P0 = 1*eye(na);  % Kalman - P(0)
Q = 10*eye(na);  % Kalman - process noise var. 
R = 1;  % Kalman - meas. noise var.

% Initialize EKF
label = 'EKF_pend';
x0 = zeros(na,1);
u0 = 0;
y0 = h(x0, u0, params);
obs = EKF_observer(na,f,h,{params},u_meas,y_meas,dfdx,dhdx,Ts,P0, ...
    Q,R,label,x0,y0);

% Check attributes
assert(obs.n == na)
assert(isequal(obs.f, @pend_StateFcn))
assert(isequal(obs.h, @pend_MeasurementFcn))
assert(isequal(obs.u_meas, u_meas))
assert(isequal(obs.y_meas, y_meas))
assert(isequal(obs.dfdx, @pend_F))
assert(isequal(obs.dhdx, @pend_H))
assert(obs.Ts == Ts)
assert(isequal(obs.params, {params}))
assert(isequal(obs.P0, P0))
assert(isequal(obs.Q, Q))
assert(isequal(obs.R, R))
assert(isequal(obs.label, label))
assert(isequal(obs.xkp1_est, x0))
assert(isequal(obs.ykp1_est, y0))

% Time series
k = (0:N-1)';
t = Ts*k;

% Initial system state
xk = bench_sim_results{1, {'x1(k+1)', 'x2(k+1)'}}';

% matrix to save the data (u and y)
data = nan(N,14);
for j = 1:N

    % Pendulum output measurement
    yk = bench_sim_results{j, 'y(k)'}';

    % Control input
    uk = bench_sim_results{j, 'u(k)'}';

    % Unmeasured input disturbance
    pk = bench_sim_results{j, 'p(k)'}';

    obs = update_EKF(obs, yk, uk);

%     % Extended Kalman filter
%     % Linearize at current operating point by
%     % re-computing Jacobians F(k) and H(k)
%     [F, H] = pend_jacob(xef,T,K,m,L,g);
%     % Compute observer gains and update covariance matrix
%     [Kf, P] = ekf_update(P,F,H,Q,R);
% 
%     % Compute next state estimates
%     % (uses non-linear model)
%     xef = [pend_xkp1(u(1),xef(1:2),K,m,L,g,T); xef(3)] + Kf*(y - H*xef);

    xef = obs.xkp1_est;
    Kf = obs.K;
    trP = trace(obs.P);

    % Record data
    data(j,:) = [k(j) t(j) uk pk yk xk' xef' Kf' trP];

    % Get states at next sample time
    xk = bench_sim_results{j, {'x1(k+1)', 'x2(k+1)'}}';

end

col_names = {'k', 't', 'u(k)', 'p(k)', 'y(k)', 'x1(k)', 'x2(k)', ...
    'xe1(k)', 'xe2(k)', 'xe3(k)', 'Kf1(k)', 'Kf2(k)', 'Kf3(k)', 'trP'};
sim_results = array2table(data, 'VariableNames', col_names);

% % Show selected results
% j = find(t == 4.5);
% selection = (0:9)' + j;
% sim_results(selection, :)
% 
% % Compare gains with benchmark data
% [
%     sim_results(selection, {'t', 'Kf1(k)', 'Kf2(k)', 'Kf3(k)', 'trP'}) ...
%     bench_sim_results(selection, {'Kf1_k_', 'Kf2_k_', 'Kf3_k_', 'trP_k_'})
% ]
% 
% % Compare estimates with benchmark data
% [
%     sim_results(selection, {'t', 'xe1(k)', 'xe2(k)', 'xe3(k)'}) ...
%     bench_sim_results(selection, {'xef1_k_1_', 'xef2_k_1_', 'xef3_k_1_'})
% ]

% Note: Homework uses a slightly different EKF calculation so estimates
% are close but not very close.
assert(isequal( ...
    round(sim_results{:, {'t', 'xe1(k)', 'xe2(k)', 'xe3(k)'}}, 1), ...
    round(bench_sim_results{:, {'t', 'xef1(k+1)', 'xef2(k+1)', 'xef3(k+1)'}}, 1) ...
))
% Covariance matrices are very different
% assert(isequal( ...
%     round(sim_results{:, {'t', 'trP'}}, 1), ...
%     round(bench_sim_results{:, {'t', 'trP(k)'}}, 1) ...
% ))


%% Test with arom3 process model

addpath('~/process-models/arom3/')

% Load non-linear model and parameters
arom3_params

% System dimensions
n = 3;
ny = 2;
nu = 2;
assert(size(x0, 1) == n)
assert(size(p0, 1) == nu)

% Augmented system dimensions
na = n + 2;

% Time step for observer
Ts = 1/60;  % hours

% Observer parameters
f = @arom3_StateFcnRodin;
h = @arom3_MeasurementFcnRodin2;
u_meas = [false; false];
y_meas = [true; true];
dfdx = @arom3_StateJacobianFcnRodin;
dhdx = @arom3_MeasurementJacobianFcnRodin2;
P0 = diag([1 25 25 0.5 0.5]);  % see p268
R = diag([1; 100]);  % see p268
Q = diag([0.25; 1; 1; 0.0033; 0.0033]);
label = 'EKF1';
x0 = [742; 463; 537; 5; 6];
y0 = x0([1 3]);
obs = EKF_observer(na,f,h,{params},u_meas,y_meas,dfdx,dhdx,Ts,P0, ...
    Q,R,label,x0,y0);

% Check attributes
assert(obs.n == na)
assert(isequal(obs.f, @arom3_StateFcnRodin))
assert(isequal(obs.h, @arom3_MeasurementFcnRodin2))
assert(isequal(obs.u_meas, u_meas))
assert(isequal(obs.y_meas, y_meas))
assert(isequal(obs.dfdx, @arom3_StateJacobianFcnRodin))
assert(isequal(obs.dhdx, @arom3_MeasurementJacobianFcnRodin2))
assert(obs.Ts == Ts)
assert(isequal(obs.params, {params}))
assert(isequal(obs.P0, P0))
assert(isequal(obs.Q, Q))
assert(isequal(obs.R, R))
assert(isequal(obs.label, label))
assert(isequal(obs.xkp1_est, x0))
assert(isequal(obs.ykp1_est, y0))

% Test instantiation with unspecified initial conditions
obs_0 = EKF_observer(na,f,h,{params},u_meas,y_meas,dfdx,dhdx,Ts,P0, ...
    Q,R,label,x0);
assert(isequal(obs_0.ykp1_est, zeros(2, 1)))
obs_0 = EKF_observer(na,f,h,{params},u_meas,y_meas,dfdx,dhdx,Ts,P0, ...
    Q,R,label);
assert(isequal(obs_0.xkp1_est, zeros(5, 1)))
assert(isequal(obs_0.ykp1_est, zeros(2, 1)))

% Load simulation data for testing observers
filename = 'arom3_sim_benchmark.csv';
data_dir = 'results';
benchmark_sim_data = readtable(fullfile(data_dir,filename), ...
    'PreserveVariableNames',true);
Y_m = benchmark_sim_data{:, {'y1_m', 'y2_m'}};

nT = size(Y_m, 1) - 1;

% Do observer calculations
sim_data = nan(nT+1, 2+na);
k_ind = (0:nT)';
t = Ts*k_ind;
for i = 1:numel(k_ind)
    k = k_ind(i);

    % First measurement y_m(0)
    yk_m = Y_m(i,:)';

    % This system has no manipulatable inputs
    uk = [];

    % Update observer
    obs = update_EKF(obs, yk_m, uk, Ts);

    % Save results
    sim_data(i,:) = [k t(i) obs.xkp1_est'];

end

x_est_labels = {'xpred1', 'xpred2', 'xpred3', 'xpred4', 'xpred5'};
assert(all(abs(sim_data(:, 3:7) - benchmark_sim_data{:, x_est_labels}) < 1e-10, [1 2]))
