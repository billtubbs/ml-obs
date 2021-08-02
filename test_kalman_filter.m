% Test functions kalman_filter.m and update_KF.m

clear all

% SISO system example from GEL-7029 course, homework 12.
% See file /gel-7029/homework/hw12/hw12_p3_kalman.m

sys_test_siso2

% Check if benchmark simulation data file exists
if ~isfile('results/hw12_p3_kalman_sim_benchmark.csv')
    error("Run 'kalman_benchmark_hw12_p3.m' to generate benchmark data.")
end


%% Define and simulate Kalman filters

% Observer parameters
W = 0.5; % estimate of Wp used in the filter design
V = 0.8; % estimate of Vp used in the filter design
P0 = eye(n)*1000; % Initialize covariance matrix
Q = diag(repmat(W,n,1));
R = diag(repmat(V,ny,1));

% Define dynamic Kalman filter using kalman_filter function
KF = kalman_filter(A,B,C,D,Ts,P0,Q,R,"KF");

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

% Intialize system (at k = 0)
x0 = zeros(n, 1);
x = x0;

% Input signal - pseudo-random binary sequence
warnId = 'Controllib:estimation:initialPRBSSequence';
warnStruct = warning('off',warnId);
U = idinput(nT+1, 'PRBS', [0 0.5]);
warning(warnStruct);

% Matrices to collect simulation data
xNprocess = zeros(n, nT+1); % process states
yNprocess = zeros(ny, nT+1); % process outputs
xNkalman = zeros(n, nT+1); % estimated states
yNkalman = zeros(ny, nT+1); % estimated process outputs
KNkalman = zeros(n, nT+1); % observer correction gains
diagPNkalman = zeros(n, nT+1); % diag of observer covariance matrix

t = Ts * (0:nT)';

for i = 1:nT

    % Process output in current timestep
    y = C*x + v(i);

    % Record process states and output
    xNprocess(:, i) = x;
    yNprocess(:, i) = y;

    % Record Kalman filter estimates of current
    % states and process outputs (made in previous
    % timestep)
    xNkalman(:, i) = KF.xkp1_est;
    yNkalman(:, i) = C * KF.xkp1_est;

    % Update KF
    KF = update_KF(KF, U(i), y);

    % Record Kalman filter variables
    KNkalman(:, i) = KF.K;
    diagPNkalman(:, i) = diag(KF.P);

    % Process states in next timestep
    x = A*x + B*U(i) + w(:,i);

end

% Record final Kalman filter estimates
xNkalman(:, nT) = KF.xkp1_est;
yNkalman(:, nT) = C * KF.xkp1_est;


sim_results = [table(t,U) ...
    array2table(KNkalman', 'VariableNames', {'K1', 'K2'}) ...
    array2table(diagPNkalman', 'VariableNames', {'P1', 'P2'}) ...
    array2table(xNprocess', 'VariableNames', {'x1', 'x2'}) ...
    array2table(xNkalman', 'VariableNames', {'x1_est_KF', 'x2_est_KF'}) ...
    array2table(yNprocess', 'VariableNames', {'y'}) ...
    array2table(yNkalman', 'VariableNames', {'y_est_KF'})];

% Display results
%head(sim_results)


% Verify results by comparing with Kalman_Filter.mlx

filename = 'hw12_p3_kalman_sim_benchmark.csv';

warnId = 'MATLAB:table:ModifiedAndSavedVarnames';
warnStruct = warning('off',warnId);
bench_sim_results = readtable(fullfile('results', filename));
warning(warnStruct);

%head(bench_sim_results)

% Check state estimates
assert(isequal( ...
    round(sim_results{1:10, {'x1', 'x2'}}, 6), ...
    round(bench_sim_results{1:10, {'X_t__1', 'X_t__2'}}, 6) ...
))

% Check correction gains
assert(isequal( ...
    round(sim_results{1:10, {'K1', 'K2'}}, 6), ...
    round(bench_sim_results{1:10, {'K_t__1', 'K_t__2'}}, 6) ...
))

% Check P covariances
assert(isequal( ...
    round(sim_results{1:10, {'P1', 'P2'}}, 6), ...
    round(bench_sim_results{1:10, {'P_t__1', 'P_t__2'}}, 6) ...
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

