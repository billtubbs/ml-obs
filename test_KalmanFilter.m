% Test KalmanFilter and KalmanFilterSS classes
% (See class definitions in KalmanFilter.m and KalmanFilterSS.m)

clear all

% Folder containing test data
results_dir = 'results';


%% SISO system example from GEL-7029
% See file Kalman_Filter.mlx

sys_test_siso

% Check if benchmark simulation data file exists
filename = 'KF_sim_benchmark.csv';
if ~isfile(fullfile(results_dir, filename))
    error("Test data file '%s' not found.\n Run 'Kalman_Filter_benchmark.mlx' to generate benchmark data.", filename)
end
bench_sim_results = readtable(fullfile(results_dir, filename));

% Covariance matrices
Q = diag([0.1; 0.1]);
R = 0.5;
N = zeros(n,ny);
x0 = [0.1; 0.5];

% Define steady-state Kalman filter with KalmanFilterSS class
label = 'KFSS';
KFSS_old = kalman_filter_ss(A,B,C,D,Ts,Q,R,label,x0);
KFSS = KalmanFilterSS(A,B,C,D,Ts,Q,R,label,x0);

assert(strcmp(KFSS.type, "KFSS"))
assert(isequal(KFSS.A, A))
assert(isequal(KFSS.B, B))
assert(isequal(KFSS.C, C))
assert(isequal(KFSS.D, D))
assert(isequal(KFSS.Ts, Ts))
assert(isequal(KFSS.Q, Q))
assert(isequal(KFSS.R, R))
assert(isequal(KFSS.label, label))
assert(isequal(KFSS.xkp1_est, x0))
assert(KFSS.ykp1_est == C * x0)
assert(KFSS.n == n)
assert(KFSS.nu == nu)
assert(KFSS.ny == ny)

% Re-define with no initial state specified (should be set to zero)
KFSS = KalmanFilterSS(A,B,C,D,Ts,Q,R,"KFSS");
K_test = [0.7727; 0.7557];
P_test = [1.5098    1.2170;
          1.2170    1.2191];

assert(isequal(round(KFSS.K, 4), K_test))
assert(isequal(round(KFSS.P, 4), P_test))
assert(isequal(KFSS.xkp1_est, zeros(n, 1)))
assert(KFSS.ykp1_est == 0)

% Define dynamic Kalman filter with KalmanFilter class
P0 = diag([1e-4 1e-4]);
KF = KalmanFilter(A,B,C,D,Ts,P0,Q,R,"KF");

assert(all(isnan(KF.K)))
assert(isequal(KF.P, P0))
assert(isequal(KF.xkp1_est, zeros(2, 1)))
assert(KF.ykp1_est == 0)

% number of points to simulate
nT = 100;

% seed random number generator
rng(0)

% Measurement noise for the whole simulation
v = sqrt(Rp)*randn(nT,1);

% Process noise for the whole simulation
w = sqrt(Qp)*randn(2,nT); 

u0 = 1;  % initial value of u
x0 = (eye(length(A)) - A) \ (B*u0);  % steady-state value of x

% Intialize system (at k = 0)
x = x0;

% Input signal
U = [zeros(10,1); ones(nT+1-10, 1)]; % process input for the whole simulation

% Matrices to collect simulation data
xNprocess = zeros(n, nT+1); % process states
yNprocess = zeros(ny, nT+1); % process outputs
xNkalman1 = zeros(n, nT+1); % estimated states
yNkalman1 = zeros(ny, nT+1); % estimated process outputs
xNkalman2 = zeros(n, nT+1); % estimated states
yNkalman2 = zeros(ny, nT+1); % estimated process outputs

for i = 1:nT

    % Process output in current timestep
    y = C*x + v(i);
    
    % Record process states and output
    xNprocess(:, i) = x;
    yNprocess(:, i) = y; 

    % Process states in next timestep
    x = A*x + B*U(i) + w(:,i);

    % Update Kalman filters
    KFSS.update(U(i), y);
    KF.update(U(i), y);

    % Record Kalman filter estimates in next timestep
    xNkalman1(:, i+1) = KF.xkp1_est;
    yNkalman1(:, i+1) = KF.ykp1_est;
    xNkalman2(:, i+1) = KFSS.xkp1_est;
    yNkalman2(:, i+1) = KFSS.ykp1_est;

end
t = Ts * (0:nT)';


% plot results

% figure(1); clf
% 
% subplot(411);
% plot(t', yNprocess, 'k', t', yNkalman1, 'r', t', yNkalman2, 'g', 'Linewidth', 2)
% legend('Process output', 'KF1 estimates', 'KF2 estimates')
% ylabel('y_1')
% grid on
% 
% subplot(412);
% plot(t', xNprocess(1,:), 'k', t', xNkalman1(1,:), 'r', ...
%     t', xNkalman2(1,:), 'g', 'Linewidth', 2)
% legend('Process state', 'KF1 estimates', 'KF2 estimates')
% ylabel('x_1')
% grid on
% 
% subplot(413);
% plot(t', xNprocess(2,:), 'k', t', xNkalman1(2,:), 'r', ...
%     t', xNkalman2(2,:), 'g', 'Linewidth', 2);
% legend('Process state', 'KF1 estimates', 'KF2 estimates')
% ylabel('x_2')
% grid on
% 
% subplot(414);
% stairs(t', U', 'Linewidth', 2);
% xlabel('Time [s]');
% ylabel('u_1')
% grid on

% Display results
sim_results = [table(t,U) ...
    array2table(xNprocess', 'VariableNames', {'x1', 'x2'}) ...
    array2table(xNkalman1', 'VariableNames', {'x1_est_KF', 'x2_est_KF'}) ...
    array2table(xNkalman2', 'VariableNames', {'x1_est_KFSS', 'x2_est_KFSS'}) ...
    array2table(yNprocess', 'VariableNames', {'y'}) ...
    array2table(yNkalman1', 'VariableNames', {'y_est_KF'}) ...
    array2table(yNkalman2', 'VariableNames', {'y_est_KFSS'})];

%head(sim_results)

% Verify results by comparing with outputs of Kalman_Filter.mlx

%head(bench_sim_results)

assert(isequal( ...
    round(sim_results{1:100, {'x1', 'x2'}}, 7), ...
    round(bench_sim_results{1:100, {'xNprocess_1', 'xNprocess_2'}}, 7) ...
))


%% Test on 2x2 system

% Check if benchmark simulation data file exists
filename = 'KF_sim_benchmark_2x2.csv';
if ~isfile(fullfile(results_dir, filename))
    error("Test data file '%s' not found.", filename)
end

bench_sim_results = readtable(fullfile(results_dir, filename));

% Set RNG seed
rng(0)

% Noise variances
sigma_p = 0;  % input disturbances
sigma_M = 0.1;  % measurement noise

% Sample time
Ts = 1;

% Discrete time state space model
A = [ 0.8890       0     1 -0.2;
           0  0.8890  -0.2    1;
           0       0     1    0;
           0       0     0    1];
B = [    1 -0.2  0  0;  % TODO: increase the coupling, -0.5?
      -0.2    1  0  0;
         0    0  1  0;
         0    0  0  1];
C = [ 0.1110 0         0  0;
             0  0.1110 0  0];
D = zeros(2, 4);
Gpss = ss(A,B,C,D,Ts);

% Dimensions
n = size(A, 1);
nu = size(B, 2);
ny = size(C, 1);

% Designate known inputs and measured outputs
u_meas = [1 2];
y_meas = [1 2];
np = nu - numel(u_meas);
nu = nu - np;
Bu = B(:, u_meas);
Du = D(:, u_meas);

% Default initial condition
x0 = zeros(n, 1);

% Observer parameters
Q = 0.01 .* eye(n);
R = 0.1 .* eye(ny);
P0 = 1e-4 .* eye(n);

% Define Kalman filters
KFSS_test = kalman_filter_ss(A,Bu,C,Du,Ts,Q,R,"KFSS",x0);
KFSS = KalmanFilterSS(A,Bu,C,Du,Ts,Q,R,"KFSS",x0);
KF_test = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,"KF",x0);
KF = KalmanFilter(A,Bu,C,Du,Ts,P0,Q,R,"KF",x0);
n_obs = 4;

% Number of timesteps to simulate
nT = 100;
t = Ts*(0:nT)';

% Random inputs
P = sigma_p^2 .* [zeros(5, np); randn([45 np]); zeros(nT+1-50, np)];
U = [zeros(5, nu); idinput([45 nu]); [-0.1 0.25].*ones(nT+1-50, nu)];
V = sigma_M^2 .* randn([nT+1 ny]);

% Simulate system
U_sim = [U P];
[Y, t, X] = lsim(Gpss, U_sim, t, x0);

% Add measurement noise
Y_m = Y + V;

% Plot results

% figure(2); clf
% 
% ax1 = subplot(311);
% plot(t, Y_m(:, 1), 'o', t, Y_m(:, 2), 'o'); hold on
% % Modify plot colors
% cols = get(gca,'ColorOrder');
% cols(ny+1:2*ny, :) = cols(1:ny, :);
% set(gca,'ColorOrder', cols);
% plot(t, Y(:, 1), t, Y(:, 2), 'Linewidth', 2)
% legend({'$y_{m,1}(k)$', '$y_{m,2}(k)$', '$y_1(k)$', '$y_2(k)$'}, 'Interpreter', 'latex')
% ylabel('$y_i(k)$', 'Interpreter', 'latex')
% title('Outputs')
% grid on
% 
% ax2 = subplot(312);
% for i = 1:nu
%     stairs(t, U(:, i), 'Linewidth', 2); hold on
% end
% ylim([-1.2 1.2])
% legend({'$u_1(k)$', '$u_2(k)$'}, 'Interpreter', 'latex')
% ylabel('$u_i(k)$', 'Interpreter', 'latex')
% title('Known inputs')
% grid on
% 
% ax3 = subplot(313);
% for i = 1:np
%     stairs(t, P(:, i), 'Linewidth', 2); hold on
% end
% legend({'$p_1(k)$', '$p_2(k)$'}, 'Interpreter', 'latex')
% ylabel('$p_i(k)$', 'Interpreter', 'latex')
% title('Unknown inputs')
% grid on

Xkp1_est = zeros(nT+1, n*n_obs);
Ykp1_est = zeros(nT+1, ny*n_obs);
E_obs = zeros(nT+1, ny*n_obs);
trP_obs = zeros(nT+1, n_obs);

for i = 1:nT

    yk = Y_m(i, :)';
    if i > 1
        uk = U(i-1, :)';
    else
        uk = zeros(nu, 1);
    end

    % Update observer gains and covariance matrix
    KFSS.update(yk, uk);
    KF.update(yk, uk);
    KFSS_test = update_KF(KFSS_test, uk, yk);
    KF_test = update_KF(KF_test, uk, yk);

    % Record observer 1 estimates
    Xkp1_est(i, 1:n) = KFSS.xkp1_est';
    Ykp1_est(i, 1:ny) = KFSS.ykp1_est';
    E_obs(i, 1:ny) = yk' - KFSS.ykp1_est';
    trP_obs(i, 1) = trace(KFSS.P);

    % Record observer 2 estimates
    Xkp1_est(i, n+1:2*n) = KF.xkp1_est';
    Ykp1_est(i, ny+1:2*ny) = KF.ykp1_est';
    E_obs(i, ny+1:2*ny) = yk' - KF.ykp1_est';
    trP_obs(i, 2) = trace(KF.P);

    % Record observer 3 estimates
    Xkp1_est(i, 2*n+1:3*n) = KFSS_test.xkp1_est';
    Ykp1_est(i, 2*ny+1:3*ny) = KFSS_test.ykp1_est';
    E_obs(i, 2*ny+1:3*ny) = yk' - KFSS_test.ykp1_est';
    trP_obs(i, 3) = trace(KFSS_test.P);

    % Record observer 4 estimates
    Xkp1_est(i, 3*n+1:4*n) = KF_test.xkp1_est';
    Ykp1_est(i, 3*ny+1:4*ny) = KF_test.ykp1_est';
    E_obs(i, 3*ny+1:4*ny) = yk' - KF_test.ykp1_est';
    trP_obs(i, 4) = trace(KF_test.P);

end

% Save benchmark results - from observers 3, 4
% bench_sim_results = [table(t,U,P,Y,V,Y_m) ...
%     array2table(X, 'VariableNames', {'X_1', 'X_2', 'X_3', 'X_4'}) ...
%     array2table(Xkp1_est(:, 2*n+1:4*n), 'VariableNames', {'X_est_1', 'X_est_2', 'X_est_3', 'X_est_4', 'X_est_5', 'X_est_6', 'X_est_7', 'X_est_8'}) ...
%     array2table(Ykp1_est(:, 2*ny+1:4*ny), 'VariableNames', {'Y_est_1', 'Y_est_2', 'Y_est_3', 'Y_est_4'}) ...
%     array2table(trP_obs(:, 3:4), 'VariableNames', {'trP_obs_1', 'trP_obs_2'})];
% writetable(bench_sim_results, fullfile(results_dir, filename))

% Combine results into table
sim_results = [table(t,U,P,Y,V,Y_m) ...
    array2table(X, 'VariableNames', {'X_1', 'X_2', 'X_3', 'X_4'}) ...
    array2table(Xkp1_est(:, 1:2*n), 'VariableNames', {'X_est_1', 'X_est_2', 'X_est_3', 'X_est_4', 'X_est_5', 'X_est_6', 'X_est_7', 'X_est_8'}) ...
    array2table(Ykp1_est(:, 1:2*ny), 'VariableNames', {'Y_est_1', 'Y_est_2', 'Y_est_3', 'Y_est_4'}) ...
    array2table(trP_obs(:, 1:2), 'VariableNames', {'trP_obs_1', 'trP_obs_2'})];
%sim_results

% Plot observer estimates
% figure(3); clf
% plot(t, Ykp1_est(:, 1), 'o', t, Ykp1_est(:, 2), 'o'); hold on
% % Modify plot colors
% cols = get(gca,'ColorOrder');
% cols(ny+1:2*ny, :) = cols(1:ny, :);
% set(gca,'ColorOrder', cols);
% plot(t, Y(:, 1), t, Y(:, 2), 'Linewidth', 2)
% xlim(t([1 end]))
% legend({'$\hat{y}_1(k)$', '$\hat{y}_2(k)$', '$y_1(k)$', '$y_2(k)$'}, 'Interpreter', 'latex')
% ylabel('$y_i(k)$', 'Interpreter', 'latex')
% title('Outputs')
% grid on

% Verify results by comparing with saved benchmark results
%head(bench_sim_results)
assert(isequal( ...
    round(sim_results.Variables, 7), ...
    round(bench_sim_results.Variables, 7) ...
))
