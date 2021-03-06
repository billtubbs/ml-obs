% Test functions MEKF_observer.m and update_MEKF.m
%
% To run this test use the following command:
%
% >> runtests MEKF_observer
%

clear all

addpath('~/process-models/pend/')

seed = 1;
rng(seed)

% Sample period
Ts = 0.1;

% Noise variances
sigma_W = [0; 0; 0];
sigma_M = 0.01;

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
state_fcn = @pend_StateFcn;
meas_fcn = @pend_MeasurementFcn;
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


% Define observers

% Define extended Kalman filter 1
x0 = [0.1; -0.5; 0];
u0 = 0;
y0 = meas_fcn(x0,u0,params1);
EKF1 = EKF_observer(na,state_fcn,meas_fcn,{params1},u_meas,y_meas,dfdx,dhdx,Ts,P0, ...
    Q1,R1,'EKF1',x0,y0);
assert(isequal(EKF1.xkp1_est, x0))
assert(isequal(EKF1.ykp1_est, y0))

% Define extended Kalman filter 2
x0 = [0.1; -0.5; 0];
u0 = 0;
y0 = meas_fcn(x0,u0,params1);
EKF2 = EKF_observer(na,state_fcn,meas_fcn,{params2},u_meas,y_meas,dfdx,dhdx,Ts,P0, ...
    Q2,R2,'EKF1',x0,y0);
assert(isequal(EKF2.xkp1_est, x0))
assert(isequal(EKF2.ykp1_est, y0))

% Transition probabilities
T = [0.9 0.1; 0.1 0.9];
assert(all(sum(T, 2) == 1))

% System indicator sequences
seq = {
    [0 0 0 0 0]; ...
    [0 0 0 0 1]; ...
    [0 0 0 1 1]; ...
    [0 0 1 1 1]; ...
    [0 1 1 1 1]; ...
    [1 1 1 1 1]; ...
    [1 1 1 1 0]; ...
    [1 1 1 0 0]; ...
    [1 1 0 0 0]; ...
    [1 0 0 0 0] ...
    };
n_filt = numel(seq);

% Define multi-model Kalman filter
state_fcn = {state_fcn, state_fcn};
meas_fcn = {meas_fcn, meas_fcn};
params = {{params1}, {params2}};
dfdx = {dfdx, dfdx};
dhdx = {dhdx, dhdx};
Q = {Q1, Q2};
R = {R1, R2};
P0 = repmat({P0}, n_filt, 1);
x0 = [0.1; -0.5; 0];  % TODO: Check x0 is working
u0 = 0;
y0 = meas_fcn{1}(x0,u0,params1);
d = 1;
MEKF = MEKF_observer(na,state_fcn,meas_fcn,params,u_meas,y_meas,dfdx,dhdx,Ts, ...
    P0,Q,R,seq,T,d,"MEKF",x0,y0);
assert(isequal(MEKF.xkp1_est, x0))
assert(isequal(MEKF.ykp1_est, y0))

% Re-define with no initial state specified (should be set to zero)
MEKF = MEKF_observer(na,state_fcn,meas_fcn,params,u_meas,y_meas,dfdx,dhdx,Ts, ...
    P0,Q,R,seq,T,d,"MEKF");

assert(isequal(MEKF.state_fcn, state_fcn))
assert(isequal(MEKF.meas_fcn, meas_fcn))
assert(isequal(MEKF.params, params))
assert(isequal(MEKF.dfdx, dfdx))
assert(isequal(MEKF.dhdx, dhdx))
assert(MEKF.Ts == Ts)
assert(isequal(MEKF.P0, P0))
assert(isequal(MEKF.Q, Q))
assert(isequal(MEKF.R, R))
assert(MEKF.n_filt == n_filt)
assert(MEKF.n_filt == numel(MEKF.filters))
assert(isequaln(MEKF.i, nan(1, 2)))
assert(MEKF.n == na)
assert(MEKF.nu == nu)
assert(MEKF.ny == ny)
assert(MEKF.Ts == Ts)
assert(MEKF.f == size(MEKF.seq{1}, 2))
assert(MEKF.nj == 2)
assert(isequal(MEKF.T, T))
for j = 1:n_filt
    assert(isequal(MEKF.filters{j}.xkp1_est, zeros(na, 1)))
    assert(MEKF.filters{j}.ykp1_est == 0)
end
assert(isequal(MEKF.xkp1_est, zeros(na, 1)))
assert(MEKF.ykp1_est == 0)

% Run simulation

filename = 'pend_sim_benchmark.csv';
results_dir = 'results';
bench_sim_results = readtable(fullfile(results_dir, filename), ...
    'PreserveVariableNames',true);
% In MATLAB 2020b 'PreserveVariableNames' is deprecated.  Use
% readtable(filename,'PreserveVariableNames',true)

% Simulation parameters
N = 150;  % no. of simulation points = 15 sec
Ts = 0.1;  % sampling period

% Time series
k = (0:N-1)';
t = Ts*k;

% Input signal - base initially on benchmark simulation
U = bench_sim_results{:, 'u(k)'};
U(t>=8) = idinput(sum(t >= 8));

% Measurement noise
V = sigma_M*randn(N, ny);

% System model selection indicator
M = ones(N, 1);
M(t>=10) = 2;

% Initial system state
xk = zeros(2, 1);
col_names = {'k', 't', 'u(k)', 'p(k)', 'm(k)', 'y(k)', ...
    'x1(k)', 'x2(k)', ...
    'xe1(k) EKF1', 'xe2(k) EKF1', 'xe3(k) EKF1', ...
    'xe1(k) EKF2', 'xe2(k) EKF2', 'xe3(k) EKF2', ...
    'xe1(k) MKF', 'xe2(k) MKF', 'xe3(k) MKF'};
%fprintf("%4s  %4s  %9s  %9s  %9s  %9s  %9s  %9s  %9s  %9s\n", col_names{:});

% matrix to save the data (u and y)
data = nan(N,17);
for j = 1:N

    % Control input
    uk = U(j, :)';

    % Disturbance
    pk = bench_sim_results{j, 'p(k)'}';

    % Pendulum output measurement
    mk = M(j, 1);
    yk = pend_yk(xk, uk+pk, params{mk}{:});

    yk_m = yk + V(j);

    % Update observers
    EKF1 = update_EKF(EKF1, yk_m, uk);
    EKF2 = update_EKF(EKF2, yk_m, uk);
    MEKF = update_MEKF(MEKF, yk_m, uk);
    xef = MEKF.xkp1_est;
    xef1 = EKF1.xkp1_est;
    xef2 = EKF2.xkp1_est;
    %trP = trace(obs.P);  % TODO: Calculate this for MKF observers

    % Display results
    %fprintf("%4d  %4.1f  %9.4f  %9.4f  %9.4f  %9.4f  %9.4f  %9.4f  %9.4f  %9.4f\n", ...
    %    k(j), t(j), uk, pk, yk, xk', xef');

    % Record data
    data(j,:) = [k(j) t(j) uk pk mk yk xk' xef1' xef2' xef'];

    % Update system model states
    xk = pend_xkp1(xk, uk+pk, params{mk}{:});

end

sim_results = array2table(data, 'VariableNames', col_names);

% % Show selected results
% j = find(t == 4.5);
% selection = (0:9)' + j;
% sim_results(selection, :)

%TODO: Need to test the results somehow...

% % Compare estimates with benchmark data
% [
%     sim_results(selection, {'t', 'xe1(k) MKF', 'xe2(k) MKF', 'xe3(k) MKF'}) ...
%     bench_sim_results(selection, {'xef1(k+1)', 'xef2(k+1)', 'xef3(k+1)'})
% ]
% 
% % Plot state estimates
% figure(1); clf
% x_labels = {'angle', 'angular velocity', 'input disturbance model state'};
% x_labels = {'x1(k)', 'x2(k)', 'p(k)'};
% axes = nan(1, na);
% for i = 1:na
%     axes(i) = subplot(na, 1, i);
%     label1 = sprintf('xe%d(k) EKF1', i);
%     plot(t, sim_results{:, {label1}}, 'Linewidth', 2); hold on
%     label2 = sprintf('xe%d(k) EKF2', i);
%     plot(t, sim_results{:, {label2}}, 'Linewidth', 2);
%     label3 = sprintf('xe%d(k) MKF', i);
%     plot(t, sim_results{:, {label3}}, '--', 'Linewidth', 2);
%     label4 = x_labels{i};
%     plot(t, sim_results{:, {label4}}, 'k--');
%     ylabel(x_labels{i}, 'Interpreter','Latex');
%     grid on
%     title(sprintf('State $x_%d$', i), 'Interpreter','Latex')
%     if i == na
%         xlabel('t (hours)', 'Interpreter','Latex');
%     end
%     legend({'EKF1', 'EKF2', 'MEKF', 'true'})
% end
% 
% linkaxes(axes, 'x')

% % Plot trace of covariance matrices
% figure(2); clf
% plot(t_sim, sim_results{:, {'trP_EKF', 'trP_EKF2', 'trP_EKF3'}},'Linewidth',2)
% xlabel('t (hours)', 'Interpreter','Latex');
% ylabel('tr $P(k)$', 'Interpreter','Latex');
% set(gca, 'YScale', 'log')
% grid on
% title('Trace of Covariance Matrix', 'Interpreter','Latex')
% legend(obs_labels)