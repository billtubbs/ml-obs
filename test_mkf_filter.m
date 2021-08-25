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
MKF = mkf_filter(A,B,C,D,Ts,P0,Q,R,seq,T,"MKF");

assert(MKF.n_filt == n_filt)
assert(MKF.n_filt == numel(MKF.filters))
assert(MKF.i == 0)
assert(MKF.n == n)
assert(MKF.nu == nu)
assert(MKF.ny == ny)
assert(MKF.Ts == Ts)
assert(MKF.nf == size(MKF.seq{1}, 2))
assert(MKF.nj == 2)
assert(isequal(size(MKF.xkp1_est), [n 1]))
assert(isequal(size(MKF.ykp1_est), [ny 1]))
assert(isequal(MKF.T, T))

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
    
    % Kalman update equations
    % Update observer gains and covariance matrix
    MKF = update_MKF(MKF, [ukm1; 0], yk);

    % Record observer estimates
    X_est(i, :) = MKF.xkp1_est';
    Y_est(i, :) = MKF.ykp1_est';
    E_obs(i, :) = yk' - MKF.ykp1_est';
    
    for j=1:MKF.n_filt
        K_obs{i, j} = MKF.filters{j}.K';
        trP_obs{i, j} = trace(MKF.filters{j}.P);
    end
end

% Combine results
sim_results = table(t,Gamma,U,Wp,X,Y,X_est,Y_est,E_obs,K_obs,trP_obs);

% Display results
%sim_results

% Compare simulation results to saved data
results_dir = 'results';
filename = 'mkf_test_sim_results.csv';
test_sim_results = readtable(fullfile(results_dir, filename));
assert(all(abs(sim_results{:, 'X_est'} ...
    - test_sim_results{:, {'X_est_1', 'X_est_2'}}) < 1e-10, [1 2]))

% Compute mean-squared error
mse = mean((Y_est - Y).^2);
assert(round(mse, 4) == 0.0188)

% % Plot of inputs and outputs
% 
% figure(1); clf
% ax1 = subplot(4,1,1);
% stairs(t,Y); hold on
% stairs(t,Y_est,'Linewidth',2);
% max_min = [min(min([Y Y_est])) max(max([Y Y_est]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('y(k) and y_est(k)')
% title('Process output measurements')
% legend('y(k)','y_est(k)','Interpreter','none')
% grid on
% 
% ax2 = subplot(4,1,2);
% stairs(t,X); hold on
% stairs(t,X_est,'Linewidth',2);
% max_min = [min(min([X X_est])) max(max([X X_est]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('x(k)')
% labels = repmat({''}, 1, n*2);
% for i=1:n
%     labels{i} = sprintf("x_%d(k)",i);
% end
% for i=1:n
%     labels{i+n} = sprintf("x_%d_est(k)",i);
% end
% legend(labels,'Interpreter','none')
% title('States')
% grid on
% 
% ax3 = subplot(4,1,3);
% stairs(t,U,'Linewidth',2); hold on;
% stairs(t,Wp,'Linewidth',2)
% max_min = [min(min([U Wp])) max(max([U Wp]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('u(k) and w_p(k)')
% legend('u(k)', 'w_p(k)')
% title('Inputs')
% grid on
% 
% ax4 = subplot(4,1,4);
% stairs(t,Gamma,'Linewidth',2)
% max_min = [min(min(Gamma)) max(max(Gamma))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('gamma(k)')
% title('Random shock sequence')
% grid on
% 
% linkaxes([ax1, ax2 ax3 ax4], 'x')