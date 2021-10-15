% Test functions mkf_filter.m and update_MKF.m

clear all

seed = 1;
rng(seed)

% Sample period
Ts = 0.5;

% Input disturbance variance
sigma_w = 0.01;

% Noise variances
sigma_W = [0; 0];
sigma_M = [0.1; 0.2];

% Discrete time state space models

% Model #1
A1 = [0.7 1;
      0 1];
B1 = [1 0;
      0 1];
C1 = [0.3 0];
D1 = zeros(1, 2);
Gpss1 = ss(A1,B1,C1,D1,Ts);
Q1 = diag([0.01 0.01]);
R1 = sigma_M(1)^2;

% Model #2
A2 = [0.9 1;
      0 1];
B2 = [1 0;
      0 1];
C2 = [0.1 0];
D2 = zeros(1, 2);
Gpss2 = ss(A2,B2,C2,D2,Ts);
Q2 = diag([0.01 0.01]);
R2 = sigma_M(2)^2;

% Dimensions
n = size(A1, 1);
nu = size(B1, 2);
ny = size(C1, 1);

assert(isequal(size(A1), size(A2)))
assert(isequal(size(B1), size(B2)))
assert(isequal(size(C1), size(C2)))
assert(isequal(size(D1), size(D2)))
assert(isequal(size(Q1), size(Q2)))
assert(isequal(size(R1), size(R2)))


%% Kalman filter with multi-model switching system

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
A = {A1, A2};
B = {B1, B2};
C = {C1, C2};
D = {D1, D2};
Q = {Q1, Q2};
R = {R1, R2};
P0 = repmat({diag([1e-4 1e-4])}, n_filt, 1);
x0 = [0.1; 0.5];
MKF1 = mkf_filter(A,B,C,D,Ts,P0,Q,R,seq,T,"MKF1",x0);
assert(isequal(MKF1.xkp1_est, x0))
assert(isequal(MKF1.ykp1_est, C{1} * x0))

% Re-define with no initial state specified (should be set to zero)
MKF1 = mkf_filter(A,B,C,D,Ts,P0,Q,R,seq,T,"MKF1");

assert(MKF1.n_filt == n_filt)
assert(MKF1.n_filt == numel(MKF1.filters))
assert(MKF1.i == 0)
assert(MKF1.n == n)
assert(MKF1.nu == nu)
assert(MKF1.ny == ny)
assert(MKF1.Ts == Ts)
assert(MKF1.nf == size(MKF1.seq{1}, 2))
assert(MKF1.nj == 2)
assert(isequal(MKF1.T, T))
assert(isequal(MKF1.xkp1_est, zeros(n, 1)))
assert(MKF1.ykp1_est == 0)

% Testing update counter
MKF1.c = 0;
MKF1.d = 1;

% Define second observer (should produce identical estimates)
% Testing update counter
seq = {
    [0 0 0]; ...
    [0 0 1]; ...
    [0 1 0]; ...
    [1 0 0] ...
    };
n_filt = numel(seq);
MKF2 = mkf_filter(A,B,C,D,Ts,P0,Q,R,seq,T,"MKF2");
MKF2.c = 0;
MKF2.d = 2;

% Simulation settings
nT = 50;
t = Ts*(0:nT)';

% Random inputs
U = (idinput(size(t)) + 1)/2;
Wp = sigma_w * randn(size(t));
U_sim = [U Wp];

% Switching sequence
Gamma = int8(rand(nT+1, 1) > T(1, 1));

% Simulate system
X = zeros(nT+1,n);
Y = zeros(nT+1,ny);
xk = zeros(n,1);

for i=1:nT+1

    % Switch system
    j = Gamma(i) + 1;

    % Inputs
    uk = U_sim(i,:)';

    % Compute y(k)
    yk = C{j}*xk + D{j}*uk;

    % Store results
    X(i, :) = xk';
    Y(i, :) = yk';

    % Compute x(k+1)
    xk = A{j}*xk + B{j}*uk;

end


% Simulate Kalman filter

X_est = zeros(nT+1,n);
Y_est = zeros(nT+1,ny);
E_obs = zeros(nT+1,ny);
K_obs = cell(nT+1,4);
trP_obs = cell(nT+1,4);

for i=1:nT+1

    yk = Y(i,:)';
    if i > 1
        ukm1 = U_sim(i-1, 1);
    else
        ukm1 = 0;
    end

    % Update observer gains and covariance matrix
    MKF1 = update_MKF2(MKF1, [ukm1; 0], yk);
    MKF2 = update_MKF2(MKF2, [ukm1; 0], yk);

    fprintf("%3d p_seq_g_Yk: [%5.3f %5.3f %5.3f %5.3f]'\n", i, MKF2.p_seq_g_Yk')

    % Record observer estimates
    X_est(i, :) = MKF1.xkp1_est';
    Y_est(i, :) = MKF1.ykp1_est';
    E_obs(i, :) = yk' - MKF1.ykp1_est';

    for j=1:MKF1.n_filt
        K_obs{i, j} = MKF1.filters{j}.K';
        trP_obs{i, j} = trace(MKF1.filters{j}.P);
    end
end

% Combine results
sim_results = table(t,Gamma,U,Wp,X,Y,X_est,Y_est,E_obs,K_obs,trP_obs);

% Display results
sim_results(:,{'t', 'Gamma', 'U', 'Wp', 'X', 'X_est', 'Y', 'Y_est'})

% Plot of inputs and outputs

figure(1); clf
ax1 = subplot(4,1,1);
stairs(t,Y); hold on
stairs(t,Y_est,'Linewidth',2);
max_min = [min(min([Y Y_est])) max(max([Y Y_est]))];
bd = max([0.1 diff(max_min)*0.1]);
ylim(max_min + [-bd bd])
xlabel('t')
ylabel('y(k) and y_est(k)')
title('Process output measurements')
legend('y(k)','y_est(k)','Interpreter','none')
grid on

ax2 = subplot(4,1,2);
stairs(t,X); hold on
stairs(t,X_est,'Linewidth',2);
max_min = [min(min([X X_est])) max(max([X X_est]))];
bd = max([0.1 diff(max_min)*0.1]);
ylim(max_min + [-bd bd])
xlabel('t')
ylabel('x(k)')
labels = repmat({''}, 1, n*2);
for i=1:n
    labels{i} = sprintf("x_%d(k)",i);
end
for i=1:n
    labels{i+n} = sprintf("x_%d_est(k)",i);
end
legend(labels,'Interpreter','none')
title('States')
grid on

ax3 = subplot(4,1,3);
stairs(t,U,'Linewidth',2); hold on;
stairs(t,Wp,'Linewidth',2)
max_min = [min(min([U Wp])) max(max([U Wp]))];
bd = max([0.1 diff(max_min)*0.1]);
ylim(max_min + [-bd bd])
xlabel('t')
ylabel('u(k) and w_p(k)')
legend('u(k)', 'w_p(k)')
title('Inputs')
grid on

ax4 = subplot(4,1,4);
stairs(t,Gamma,'Linewidth',2)
max_min = [min(min(Gamma)) max(max(Gamma))];
bd = max([0.1 diff(max_min)*0.1]);
ylim(max_min + [-bd bd])
xlabel('t')
ylabel('gamma(k)')
title('Random shock sequence')
grid on

linkaxes([ax1, ax2 ax3 ax4], 'x')

% Compare simulation results to saved data
results_dir = 'results';
filename = 'mkf_test_sim_results.csv';
test_sim_results = readtable(fullfile(results_dir, filename));
assert(all(abs(sim_results{:, 'X_est'} ...
    - test_sim_results{:, {'X_est_1', 'X_est_2'}}) < 1e-10, [1 2]))

% Compute mean-squared error
mse = mean((Y_est - Y).^2);
assert(round(mse, 4) == 0.0188)

