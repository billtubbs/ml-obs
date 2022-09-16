% Test KalmanFilter, KalmanFilterF, and KalmanFilterSS classes
% (See class definitions in KalmanFilter.m, KalmanFilterF.m,
% and KalmanFilterSS.m)

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
label = "KFSS";
KFSS_old = kalman_filter_ss(A,B,C,D,Ts,Q,R,label,x0);
KFSS = KalmanFilterSS(A,B,C,Ts,Q,R,label,x0);

assert(strcmp(KFSS.type, "KFSS"))
assert(isequal(KFSS.A, A))
assert(isequal(KFSS.B, B))
assert(isequal(KFSS.C, C))
assert(isequal(KFSS.Ts, Ts))
assert(isequal(KFSS.Q, Q))
assert(isequal(KFSS.R, R))
assert(max(abs(KFSS.K - KFSS_old.K), [], [1 2]) < 1e-12)
assert(max(abs(KFSS.Pkp1 - KFSS_old.P), [], [1 2]) < 1e-12)
K_calc = KFSS.A * KFSS.Pkp1 * KFSS.C' * (KFSS.C * KFSS.Pkp1 * KFSS.C' + KFSS.R)^-1;
assert(max(abs(KFSS.K - K_calc)) < 1e-12)
assert(isequal(round(KFSS.K, 6), [0.772750; 0.755731]))
assert(isequal(round(KFSS.Pkp1, 6), [1.509786 1.216953; 1.216953 1.219071]))
assert(isequal(KFSS.label, label))
assert(isequal(KFSS.xkp1_est, x0))
assert(KFSS.ykp1_est == C*x0)
assert(KFSS.n == n)
assert(KFSS.nu == nu)
assert(KFSS.ny == ny)

% Re-define with no initial state specified (should be set to zero)
KFSS = KalmanFilterSS(A,B,C,Ts,Q,R,"KFSS");
K_test = [0.7727; 0.7557];
P_test = [1.5098    1.2170;
          1.2170    1.2191];

assert(isequal(round(KFSS.K, 4), K_test))
assert(isequal(round(KFSS.Pkp1, 4), P_test))
assert(isequal(KFSS.xkp1_est, zeros(n, 1)))
assert(KFSS.ykp1_est == 0)

% Define dynamic Kalman filters with KalmanFilter 
% and KalmanFilterF classes
P0 = diag([1e-4 1e-4]);

KF = KalmanFilter(A,B,C,Ts,P0,Q,R,"KF");
assert(isequal(KF.A, A))
assert(isequal(KF.B, B))
assert(isequal(KF.C, C))
assert(isequal(KF.Ts, Ts))
assert(isequal(KF.Pkp1, P0))
assert(isequal(KF.Q, Q))
assert(isequal(KF.R, R))
assert(all(isnan(KF.K)))
assert(isequal(KF.xkp1_est, zeros(2, 1)))
assert(KF.ykp1_est == 0)

% New version
model.A = A;
model.B = B;
model.C = C;
model.Q = Q;
model.R = R;
model.Ts = Ts;
KFF = KalmanFilterF(model,P0,"KF");
assert(isequal(KFF.model, model))
assert(all(isnan(KFF.Kf)))
assert(isequal(KFF.xkp1_est, zeros(2, 1)))
assert(KFF.ykp1_est == 0)
assert(isequal(KFF.Pkp1, P0))
assert(all(isnan(KFF.xk_est)))
assert(all(isnan(KFF.yk_est)))
assert(all(isequaln(KFF.Pk, nan(2))))

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
xNkalman3 = zeros(n, nT+1); % estimated states
yNkalman3 = zeros(ny, nT+1); % estimated process outputs

for i = 1:nT

    % Process output in current timestep
    y = C*x + v(i);
    
    % Record process states and output
    xNprocess(:, i) = x;
    yNprocess(:, i) = y; 

    % Process states in next timestep
    x = A*x + B*U(i) + w(:,i);

    % Check predictions in next time step
    % KFF.predict();
    % assert(all(abs(A * KFF.xk_est + B * U(i) - KFF.xkp1_est) < 1e-14))

    % Update Kalman filters
    KFSS.update(y, U(i));
    KF.update(y, U(i));
    KFF.update(y, U(i));

    % Record Kalman filter estimates in next timestep
    xNkalman1(:, i+1) = KF.xkp1_est;
    yNkalman1(:, i+1) = KF.ykp1_est;
    xNkalman2(:, i+1) = KFSS.xkp1_est;
    yNkalman2(:, i+1) = KFSS.ykp1_est;
    xNkalman3(:, i+1) = KFF.xk_est;
    yNkalman3(:, i+1) = KFF.yk_est;

    % Check updated states and predictions match
    assert(all(abs(C * KFF.xk_est - KFF.yk_est) < 1e-14))
    assert(all(abs(C * KFF.xkp1_est - KFF.ykp1_est) < 1e-14))

    % Check predictions are the same
    assert(all(abs(KF.xkp1_est - KFF.xkp1_est) < 1e-14))
    assert(all(abs(KF.ykp1_est - KFF.ykp1_est) < 1e-14))

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


%% SISO system example from GEL-7029 course, homework 12.
% See file /gel-7029/homework/hw12/hw12_p3_kalman.m

sys_test_siso2

% Check if benchmark simulation data file exists
if ~isfile('results/hw12_p3_kalman_sim_benchmark.csv')
    error("Run 'kalman_benchmark_hw12_p3.m' to generate benchmark data.")
end

% Observer parameters
W = 0.5; % estimate of Wp used in the filter design
V = 0.8; % estimate of Vp used in the filter design
P0 = eye(n)*1000; % Initialize covariance matrix
Q = diag(repmat(W,n,1));
R = diag(repmat(V,ny,1));
x0 = [0.1; 0.5];

% Define dynamic Kalman filter using kalman_filter function
label = 'KF';
KF = KalmanFilter(A,B,C,Ts,P0,Q,R,label,x0);
assert(strcmp(KF.type, "KF"))
assert(isequal(KF.A, A))
assert(isequal(KF.B, B))
assert(isequal(KF.C, C))
assert(isequal(KF.Ts, Ts))
assert(isequal(KF.P0, P0))
assert(isequal(KF.Q, Q))
assert(isequal(KF.R, R))
assert(isequal(KF.Pkp1, P0))
assert(isequal(KF.label, label))
assert(isequal(KF.xkp1_est, x0))
assert(KF.ykp1_est == C*x0)
assert(KF.n == n)
assert(KF.nu == nu)
assert(KF.ny == ny)

% New version
label = 'KFF';
model.A = A;
model.B = B;
model.C = C;
model.Q = Q;
model.R = R;
model.Ts = Ts;
KFF = KalmanFilterF(model,P0,label,x0);
assert(strcmp(KFF.type, "KFF"))
assert(isequal(KFF.model, model))
assert(isequal(KFF.label, label))
assert(isequal(KFF.xkp1_est, x0))
assert(KFF.ykp1_est == C*x0)
assert(isequal(KFF.Pkp1, P0))
assert(all(isnan(KFF.xk_est)))
assert(all(isequaln(KFF.Pk, nan(2))))
assert(all(isnan(KFF.yk_est)))
assert(KFF.n == n)
assert(KFF.nu == nu)
assert(KFF.ny == ny)

% Re-define with no initial state specified (should be set to zero)
KF = KalmanFilter(A,B,C,Ts,P0,Q,R);
assert(all(isnan(KF.K)))
assert(isequal(KF.Pkp1, P0))
assert(isequal(KF.xkp1_est, zeros(n, 1)))
assert(KF.ykp1_est == 0)

% Re-define with no initial state specified (should be set to zero)
KFF = KalmanFilterF(model,P0);
assert(all(isnan(KFF.Kf)))
assert(isequal(KFF.xkp1_est, zeros(n, 1)))
assert(isequaln(KFF.Pkp1, P0))
assert(KFF.ykp1_est == 0)
assert(all(isnan(KFF.xk_est)))
assert(isequaln(KFF.Pk, nan(2)))
assert(all(isnan(KFF.yk_est)))

% number of points to simulate
nT = 100;

% seed random number generator
rng(0)

% Measurement noise for the whole simulation
V = nan(nT+1,1);
V(1:end-1,:) = sqrt(Rp)*randn(1,nT)';

% Process noise for the whole simulation
W = nan(nT+1,2);
W(1:end-1,:) = (sqrt(Qp)*randn(2,nT))';

% Intialize system (at k = 0)
x0 = zeros(n, 1);
x = x0;

% Input signal - pseudo-random binary sequence
U = nan(nT+1,1);
warning('off')
U(1:nT,1) = idinput(nT, 'PRBS', [0 0.5]);
warning('on')

% Matrices to collect simulation data
xNprocess = nan(n, nT+1); % process states
yNprocess = nan(ny, nT+1); % process outputs
xNkalman1 = nan(n, nT+1); % estimated states
yNkalman1 = nan(ny, nT+1); % estimated process outputs
KNkalman1 = nan(n, nT+1); % observer correction gains
diagPNkalman1 = nan(n, nT+1); % diag of observer covariance matrix
xNkalman2 = nan(n, nT+1); % estimated states
yNkalman2 = nan(ny, nT+1); % estimated process outputs
KNkalman2 = nan(n, nT+1); % observer correction gains
diagPNkalman2 = nan(n, nT+1); % diag of observer covariance matrix

t = Ts * (0:nT)';

for i = 1:nT

    % Process output in current timestep
    y = C*x + V(i,:)';

    % Record process states and output
    xNprocess(:, i) = x;
    yNprocess(:, i) = y;

    % Record Kalman filter estimates of current
    % states and process outputs (made in previous
    % timestep)
    xNkalman1(:, i) = KF.xkp1_est;
    yNkalman1(:, i) = KF.ykp1_est;

    % Check predictions are the same
    assert(all(abs(KF.xkp1_est - KFF.xkp1_est) < 1e-13))
    assert(all(abs(KF.ykp1_est - KFF.ykp1_est) < 1e-14))

    % Update KFs
    KF.update(y, U(i));
    KFF.update(y, U(i));

    % Record updated estimates of current states
    % and process outputs (filtering KF only)
    xNkalman2(:, i) = KFF.xk_est;
    yNkalman2(:, i) = KFF.yk_est;

    % Record Kalman filter variables
    KNkalman1(:, i) = KF.K;
    diagPNkalman1(:, i) = diag(KF.Pkp1);
    KNkalman2(:, i) = KFF.Kf;
    diagPNkalman2(:, i) = diag(KFF.Pk);

    % Process states in next timestep
    x = A*x + B*U(i) + W(i,:)';

end

% Record final Kalman filter estimates
xNkalman1(:, nT) = KF.xkp1_est;
yNkalman1(:, nT) = C * KF.xkp1_est;

sim_results = [table(t,U,V,W) ...
    array2table(xNprocess', 'VariableNames', {'x1', 'x2'}) ...
    array2table(xNkalman1', 'VariableNames', {'x1_est_KF', 'x2_est_KF'}) ...
    array2table(yNprocess', 'VariableNames', {'y'}) ...
    array2table(yNkalman1', 'VariableNames', {'y_est_KF'}) ...
    array2table(KNkalman1', 'VariableNames', {'K1', 'K2'}) ...
    array2table(diagPNkalman1', 'VariableNames', {'P1', 'P2'})];

% Verify results by comparing with code from course homework
filename = 'hw12_p3_kalman_sim_benchmark.csv';

warnId = 'MATLAB:table:ModifiedAndSavedVarnames';
warnStruct = warning('off',warnId);
bench_sim_results = readtable(fullfile('results', filename));
warning(warnStruct);

% Compare results to benchmark results
% tail(sim_results)
% tail(bench_sim_results)

% Check states and outputs
assert(isequal( ...
    round(sim_results{1:10, {'x1', 'x2', 'y'}}, 6), ...
    round(bench_sim_results{1:10, {'X_t__1', 'X_t__2', 'y_t_'}}, 6) ...
))

% Check state and output estimates
assert(isequal( ...
    round(sim_results{1:100, {'x1_est_KF', 'x2_est_KF', 'y_est_KF'}}, 6), ...
    round(bench_sim_results{1:100, {'X_e_t__1', 'X_e_t__2', 'y_est_t_'}}, 6) ...
))

% Check correction gains
assert(isequal( ...
    round(sim_results{1:100, {'K1', 'K2'}}, 6), ...
    round(bench_sim_results{1:100, {'K_t__1', 'K_t__2'}}, 6) ...
))

% Check P covariances
assert(isequal( ...
    round(sim_results{1:100, {'P1', 'P2'}}, 6), ...
    round(bench_sim_results{1:100, {'P_t__1', 'P_t__2'}}, 6) ...
))


% plot results

% figure(1); clf
% 
% subplot(411);
% plot(t', yNprocess, 'k', t', yNkalman, 'r', 'Linewidth', 2)
% legend('Process output', 'KF estimates')
% ylabel('$y$')
% grid on
% 
% subplot(412);
% plot(t', xNprocess(1,:), 'k', t', xNkalman(1,:), 'r', 'Linewidth', 2)
% legend('Process state', 'KF estimates')
% ylabel('$x_1$')
% grid on
% 
% subplot(413);
% plot(t', xNprocess(2,:), 'k', t', xNkalman(2,:), 'r', 'Linewidth', 2);
% legend('Process state', 'KF estimates')
% ylabel('$x_2$')
% grid on
% 
% subplot(414);
% stairs(t', U', 'Linewidth', 2);
% xlabel('Time [s]');
% ylabel('$u_1$')
% grid on


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
KFSS_test = kalman_filter_ss(A,Bu,C,Du,Ts,Q,R,"KFSS_test",x0);
KFSS = KalmanFilterSS(A,Bu,C,Ts,Q,R,"KFSS",x0);
KF_test = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,"KF_test",x0);
KF = KalmanFilter(A,Bu,C,Ts,P0,Q,R,"KF",x0);

% New version
model.A = A;
model.B = Bu;
model.C = C;
model.Q = Q;
model.R = R;
model.Ts = Ts;
KFF = KalmanFilterF(model,P0,"KFF",x0);

observers = {KFSS, KF, KFSS_test, KF_test, KFF};
n_obs = numel(observers);

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

Xk_est = nan(nT+1, n*n_obs);
Yk_est = nan(nT+1, n*n_obs);
Xkp1_est = nan(nT+1, n*n_obs);
Ykp1_est = nan(nT+1, ny*n_obs);
E_obs = nan(nT+1, ny*n_obs);
trP_obs = nan(nT+1, n_obs);

for i = 1:nT

    yk = Y_m(i, :)';
    if i > 1
        uk = U(i-1, :)';
    else
        uk = zeros(nu, 1);
    end

    % Record observer estimates of state predictions
    for j = 1:n_obs
        Xkp1_est(i, (j-1)*n+1:j*n) = observers{j}.xkp1_est';
        Ykp1_est(i, (j-1)*ny+1:j*ny) = observers{j}.ykp1_est';
        if isstruct(observers{j})
            % Old struct-based observers
            trP_obs(i, j) = trace(observers{j}.P);
        else
            switch observers{j}.type
                case {"KF", "KFSS"}
                    trP_obs(i, j) = trace(observers{j}.Pkp1);
                case {"KFF"}
                    trP_obs(i, j) = trace(observers{j}.Pk);
            end
        end
    end

    % Update observer gains and covariance matrix
    KFSS.update(yk, uk);
    KF.update(yk, uk);
    KFSS_test = update_KF(KFSS_test, uk, yk);
    assert(strcmp(observers{3}.label, "KFSS_test"))
    observers{3} = KFSS_test;
    KF_test = update_KF(KF_test, uk, yk);
    assert(strcmp(observers{4}.label, "KF_test"))
    observers{4} = KF_test;
    KFF.update(yk, uk);

    % Record updated estimates of current state (of
    % filtering observers only)
    for j = 5:5
        Xk_est(i, (j-1)*n+1:j*n) = observers{j}.xk_est';
        Yk_est(i, (j-1)*ny+1:j*ny) = observers{j}.yk_est';
        E_obs(i, (j-1)*ny+1:j*ny) = yk' - observers{j}.yk_est';
        trP_obs(i, j) = trace(observers{j}.Pk);
    end

    % Calculate estimation errors
    if i > 1
        for j = 1:4  % use prior prediction estimates
            E_obs(i, (j-1)*ny+1:j*ny) = (yk' - Ykp1_est(i-1, (j-1)*ny+1:j*ny));
        end
    end
    for j = 5:5  % use updated posterior estimates
        E_obs(i, (j-1)*ny+1:j*ny) = (yk' - Yk_est(i, (j-1)*ny+1:j*ny));
    end

end

% Save benchmark results - from observers 3, 4
% bench_sim_results = [table(t,U,P,Y,V,Y_m) ...
%     array2table(X, 'VariableNames', {'X_1', 'X_2', 'X_3', 'X_4'}) ...
%     array2table(Xkp1_est(:, 2*n+1:4*n), 'VariableNames', {'X_est_1', 'X_est_2', 'X_est_3', 'X_est_4', 'X_est_5', 'X_est_6', 'X_est_7', 'X_est_8'}) ...
%     array2table(Ykp1_est(:, 2*ny+1:4*ny), 'VariableNames', {'Y_est_1', 'Y_est_2', 'Y_est_3', 'Y_est_4'}) ...
%     array2table(trP_obs(:, 3:4), 'VariableNames', {'trP_obs_1', 'trP_obs_2'})];
% writetable(bench_sim_results, fullfile(results_dir, filename))

% Combine results into table - only first two observers
sim_results = [table(t,U,P,Y,V,Y_m) ...
    array2table(X, 'VariableNames', {'X_1', 'X_2', 'X_3', 'X_4'}) ...
    array2table(Xkp1_est(:, 1:2*n), 'VariableNames', {'X_est_1', 'X_est_2', 'X_est_3', 'X_est_4', 'X_est_5', 'X_est_6', 'X_est_7', 'X_est_8'}) ...
    array2table(Ykp1_est(:, 1:2*ny), 'VariableNames', {'Y_est_1', 'Y_est_2', 'Y_est_3', 'Y_est_4'}) ...
    array2table(trP_obs(:, 1:2), 'VariableNames', {'trP_obs_1', 'trP_obs_2'})];
%sim_results

% Plot observer estimates
% j = 4;
% figure(3); clf
% plot(t, Ykp1_est(:, (j-1)*ny+1), 'o', t, Ykp1_est(:, (j-1)*ny+2), 'o'); hold on
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
assert(isequaln( ...
    round(sim_results.Variables, 7), ...
    round(bench_sim_results.Variables, 7) ...
))

% Compare KFF (filtering) estimates with KF (prediction)
assert(strcmp(observers{2}.label, "KF"))
assert(strcmp(observers{5}.label, "KFF"))

% Check magnitude of differences
y_est_diffs = max(abs(Xk_est(:, (5-1)*n+1:5*n) ...
    - Xkp1_est(:, (2-1)*n+1:2*n)));
assert(all(y_est_diffs < [0.4 0.4 0.07 0.07]))
assert(all(y_est_diffs > 0.001))


%% Test copy methods

sys_test_siso

% Covariance matrices
Q = diag([0.1; 0.1]);
R = 0.5;
N = zeros(n,ny);
x0 = [0.1; 0.5];
P0 = diag([1e-4 1e-4]);

% Define steady-state Kalman
KFSS = KalmanFilterSS(A,B,C,Ts,Q,R,"KFSS",x0);

% Define dynamic Kalman filter
KF = KalmanFilter(A,B,C,Ts,P0,Q,R,"KF");

% Define dynamic Kalman filter - filtering form
model.A = A;
model.B = B;
model.C = C;
model.Q = Q;
model.R = R;
model.Ts = Ts;
KFF = KalmanFilterF(model,P0,"KFF");

% Test handle copy
KFSS_hcopy = KFSS;
assert(isequaln(KFSS_hcopy, KFSS))  % same values
assert(KFSS_hcopy == KFSS)  % must be same object

KFSS.x0 = [0.2; 0.5];
assert(isequal(KFSS_hcopy.x0, [0.2; 0.5]))

KF_hcopy = KF;
assert(isequaln(KF_hcopy, KF))  % same values
assert(KF_hcopy == KF)  % must be same object

KF.label = "New name";
assert(isequal(KF_hcopy.label, "New name"))

KFF_hcopy = KFF;
assert(isequaln(KFF_hcopy, KFF))  % same values
assert(KFF_hcopy == KFF)  % must be same object

KFF.model.A(1, 1) = KFF.model.A(1, 1) + 0.1;
assert(isequal(KFF_hcopy.model.A(1, 1), KFF.model.A(1, 1)))

% Redefine steady-state Kalman
KFSS = KalmanFilterSS(A,B,C,Ts,Q,R,"KFSS",x0);

% Redefine dynamic Kalman filter
KF = KalmanFilter(A,B,C,Ts,P0,Q,R,"KF");

% Redefine dynamic Kalman filter
KFF = KalmanFilterF(model,P0,"KFF");

% Test true copy
KFSS_copy = KFSS.copy();
assert(isequaln(KFSS_copy, KFSS))  % same values
assert(KFSS_copy ~= KFSS)  % must not be same object

KFSS.x0 = [0.2; 0.5];
assert(~isequal(KFSS_copy.x0, [0.2; 0.5]))

KF_copy = KF.copy();
assert(isequaln(KF_copy, KF))  % same values
assert(KF_copy ~= KF)  % must not be same object

KF.label = "New name";
assert(~isequal(KF_copy.label, "New name"))

KFF_copy = KFF.copy();
assert(isequaln(KFF_copy, KFF))  % same values
assert(KFF_copy ~= KFF)  % must not be same object

KFF.model.A(1, 1) = KFF.model.A(1, 1) + 0.1;
assert(~isequal(KFF_copy.model.A, KFF.model.A(1, 1)))
