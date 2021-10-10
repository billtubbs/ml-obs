% Test functions MEKF_observer.m and update_MEKF.m
%
% To run this test use the following command:
%
% >> runtests MEKF_observer
%

addpath('~/process-models/pend/')

seed = 1;
rng(seed)

% Sample period
Ts = 0.1;

% Input disturbance variance
sigma_w = 0.01;

% Noise variances
sigma_W = [0; 0];
sigma_M = [0.1; 0.2];

% Non-linear pendulum system model

% Pendulum #1 parameters
params1.K = 1.2;
params1.m = 0.3;
params1.L = 0.4;
params1.g = 9.8;
params1.dt = Ts;

% Pendulum #2 parameters
params2.K = 1.2;
params2.m = 0.5;  % different
params2.L = 0.4;
params2.g = 9.8;
params2.dt = Ts;

% Eqns for augmented model with one integrator
% (these may be common to both observers)
na = 3;
f = @pend_StateFcn;
h = @pend_MeasurementFcn;
u_meas = true;
y_meas = true;
dfdx = @pend_F;
dhdx = @pend_H;
P0 = 1*eye(na);  % Kalman - P(0)

% Observer model #1 parameters
Q1 = 10*eye(na);  % Kalman - process noise var. 
R1 = 1;  % Kalman - meas. noise var

% Observer model #2 parameters
Q2 = 10*eye(na);  % Kalman - process noise var. 
R2 = 1;  % Kalman - meas. noise var

% System dimensions
n = 2;
nu = 1;
ny = 1;


%% Define multi-model EKF observers

% Transition probabilities
T = [0.9 0.1; 0.9 0.1];
assert(all(sum(T, 2) == 1))

% System indicator sequences
seq = {
    [0 0 0 0 0 0]; ...
    [0 0 0 0 0 1]; ...
    [0 0 0 1 0 0]; ...
    [0 1 0 0 0 0] ...
    };
n_filt = numel(seq);

% Define and simulate multi-model Kalman filter
f = {f, f};
h = {h, h};
params = {params1, params2};
dfdx = {dfdx, dfdx};
dhdx = {dhdx, dhdx};
Q = {Q1, Q2};
R = {R1, R2};
P0 = repmat({P0}, n_filt, 1);
x0 = [0.1; -0.5; 2];
u0 = 0;
y0 = h{1}(x0,u0,params1);
MEKF = MEKF_observer(na,f,h,params,u_meas,y_meas,dfdx,dhdx,Ts, ...
    P0,Q,R,seq,T,"MEKF",x0,y0);
assert(isequal(MEKF.xkp1_est, x0))
assert(isequal(MEKF.ykp1_est, y0))

% Re-define with no initial state specified (should be set to zero)
MEKF = MEKF_observer(na,f,h,params,u_meas,y_meas,dfdx,dhdx,Ts, ...
    P0,Q,R,seq,T,"MEKF");

assert(isequal(MEKF.f, f))
assert(isequal(MEKF.h, h))
assert(isequal(MEKF.params, params))
assert(isequal(MEKF.dfdx, dfdx))
assert(isequal(MEKF.dhdx, dhdx))
assert(MEKF.n_filt == n_filt)
assert(MEKF.n_filt == numel(MEKF.filters))
assert(MEKF.i == 0)
assert(MEKF.n == na)
assert(MEKF.nu == nu)
assert(MEKF.ny == ny)
assert(MEKF.Ts == Ts)
assert(MEKF.nf == size(MEKF.seq{1}, 2))
assert(MEKF.nj == 2)
assert(isequal(MEKF.T, T))
for j = 1:n_filt
    assert(isequal(MEKF.filters{j}.xkp1_est, zeros(na, 1)))
    assert(MEKF.filters{j}.ykp1_est == 0)
end
assert(isequal(MEKF.xkp1_est, zeros(na, 1)))
assert(MEKF.ykp1_est == 0)


% Run simulation
obs = MEKF;

filename = 'pend_sim_benchmark.csv';
results_dir = 'results';
bench_sim_results = readtable(fullfile(results_dir, filename), ...
    'PreserveVariableNames',true);

% Simulation parameters
N = 150;  % no. of simulation points = 15 sec
Ts = 0.1;  % sampling period

% Time series
k = (0:N-1)';
t = Ts*k;

% Initial system state
xk = zeros(2, 1);
col_names = {'k', 't', 'u(k)', 'p(k)', 'yk', 'x1(k)', 'x2(k)', ...
    'xe1(k)', 'xe2(k)', 'xe3(k)'};
%fprintf("%4s  %4s  %9s  %9s  %9s  %9s  %9s  %9s  %9s  %9s\n", col_names{:});

% matrix to save the data (u and y)
data = nan(N,10);
for j = 1:N

    % Control input
    uk = bench_sim_results{j, 'u(k)'}';

    % Pendulum output measurement
    yk = pend_yk(xk, uk, params{1});

    % Disturbance
    pk = bench_sim_results{j, 'p(k)'}';

%     % debugging
%     if k(j) == 5
%         disp('stop')
%     end

    obs = update_MEKF(obs, yk, uk);

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
%    Kf = obs.K;
%    trP = trace(obs.P);
% 
%     assert(isequal( ...
%         round(xef', 4), ...
%         round(bench_sim_results{j, {'xef1(k+1)', 'xef2(k+1)', 'xef3(k+1)'}}, 4) ...
%     ))
%     assert(round(trP, 4) == round(bench_sim_results{j, 'trP(k)'}, 4))

    % Display results
    %fprintf("%4d  %4.1f  %9.4f  %9.4f  %9.4f  %9.4f  %9.4f  %9.4f  %9.4f  %9.4f\n", ...
    %    k(j), t(j), uk, pk, yk, xk', xef');

    % Record data
    data(j,:) = [k(j) t(j) uk pk yk xk' xef'];

    % Get states at next sample time
    xk = bench_sim_results{j, {'x1(k+1)', 'x2(k+1)'}}';

end

sim_results = array2table(data, 'VariableNames', col_names);

% % Show selected results
% j = find(t == 4.5);
% selection = (0:9)' + j;
% sim_results(selection, :)
% 
% % Compare gains with benchmark data
% [
%     sim_results(selection, {'t', 'Kf1(k)', 'Kf2(k)', 'Kf3(k)', 'trP'}) ...
%     bench_sim_results(selection, {'Kf1(k)', 'Kf2(k)', 'Kf3(k)', 'trP(k)'})
% ]
% 
% % Compare estimates with benchmark data
% [
%     sim_results(selection, {'t', 'xe1(k)', 'xe2(k)', 'xe3(k)'}) ...
%     bench_sim_results(selection, {'xef1(k+1)', 'xef2(k+1)', 'xef3(k+1)'})
% ]