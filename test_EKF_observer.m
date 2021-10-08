% Test EKF_observer.m and update_EKF.m

%% Test with arom3 process model

addpath('~/process-models/arom3/')

% Load non-linear model and parameters
arom3_params

n = 5;
f = @arom3_StateFcnRodin;
h = @arom3_MeasurementFcnRodin2;
u_meas = [false; false];
y_meas = [true; true];
dfdx = @arom3_StateJacobianFcnRodin;
dhdx = @arom3_MeasurementJacobianFcnRodin2;
% see Robertson et al. (1998) for simulation details
Ts = 1;
P0 = diag([1 25 25 0.5 0.5]);
R = diag([1; 100]);
sigma_w = [0.0033; 0.0033];
Q = diag([0.1; 0.1; 0.1; sigma_w]);
label = 'EKF1';
x0 = [742; 463; 537; 5; 6];
y0 = x0([1 3]);
obs = EKF_observer(n,f,h,u_meas,y_meas,dfdx,dhdx,Ts,P0,Q,R, ...
    label,x0,y0);

% Check attributes set
assert(obs.n == n)
assert(isequal(obs.f, @arom3_StateFcnRodin))
assert(isequal(obs.h, @arom3_MeasurementFcnRodin2))
assert(isequal(obs.u_meas, u_meas))
assert(isequal(obs.y_meas, y_meas))
assert(isequal(obs.dfdx, @arom3_StateJacobianFcnRodin))
assert(isequal(obs.dhdx, @arom3_MeasurementJacobianFcnRodin2))
assert(obs.Ts == Ts)
assert(isequal(obs.P0, P0))
assert(isequal(obs.Q, Q))
assert(isequal(obs.R, R))
assert(isequal(obs.label, label))
assert(isequal(obs.xkp1_est, x0))
assert(isequal(obs.ykp1_est, y0))

% Test instantiation with unspecified initial conditions
obs_0 = EKF_observer(n,f,h,u_meas,y_meas,dfdx,dhdx,Ts,P0,Q,R, ...
    label,x0);
assert(isequal(obs_0.ykp1_est, zeros(2, 1)))
obs_0 = EKF_observer(n,f,h,u_meas,y_meas,dfdx,dhdx,Ts,P0,Q,R, ...
    label);
assert(isequal(obs_0.xkp1_est, zeros(5, 1)))
assert(isequal(obs_0.ykp1_est, zeros(2, 1)))

% Test update function
uk = [];
yk_m = x0([1 3]);
dt = 1;
obs = update_EKF(obs, yk_m, uk, dt, params);
%obs.xkp1_est
%obs.ykp1_est

%TODO: run a simulation test


%% Example from homework from GEL-7029 course

clear all
addpath('~/process-models/pend/')

filename = 'hw15_EKF_sim_benchmark.csv';
results_dir = 'results';
bench_sim_results = readtable(fullfile(results_dir, filename));

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
x0 = zeros(na,1);
u0 = 0;
y0 = h(x0, u0, params);
obs = EKF_observer(na,f,h,u_meas,y_meas,dfdx,dhdx,Ts,P0,Q,R, ...
    'EKF_pend',x0,y0);

assert(isequal(obs.xkp1_est, x0))
assert(isequal(obs.ykp1_est, y0))

% Time series
k = (0:N-1)';
t = Ts*k;

% Initial system state
xk = bench_sim_results{1, {'x1_k_1_', 'x2_k_1_'}}';

% matrix to save the data (u and y)
data = nan(N,14);
for j = 1:N

    % Pendulum output measurement
    yk = bench_sim_results{j, 'y_k_'}';

    % Control input
    uk = bench_sim_results{j, 'u_k_'}';

    % Disturbance
    pk = bench_sim_results{j, 'p_k_'}';

    obs = update_EKF(obs, yk, uk, params);

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

    assert(isequal( ...
        round(xef', 4), ...
        round(bench_sim_results{j, {'xef1_k_1_', 'xef2_k_1_', 'xef3_k_1_'}}, 4) ...
    ))
    assert(round(trP, 4) == round(bench_sim_results{j, 'trP_k_'}, 4))

    % Record data
    data(j,:) = [k(j) t(j) uk pk yk xk' xef' Kf' trP];

    % Get states at next sample time
    xk = bench_sim_results{j, {'x1_k_1_', 'x2_k_1_'}}';

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

assert(isequal( ...
    round(sim_results{:, {'t', 'xe1(k)', 'xe2(k)', 'xe3(k)'}}, 4), ...
    round(bench_sim_results{:, {'t', 'xef1_k_1_', 'xef2_k_1_', 'xef3_k_1_'}}, 4) ...
))