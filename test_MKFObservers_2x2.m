% Tests the multi-model observers on a 2x2 system with RODD disturbance
% 
%  - KalmanFilterF
%      Single Kalman Filter (used for comparison).
%  - SKFObserver
%      KF with switching system model.
%  - MKFObserver
%      Multiple-model KF with switching system models.
%


clear all

addpath("~/ml-plot-utils/")

seed = 0;
rng(seed)


%% Simulation test on 2x2 linear system with RODD disturbances

% Sample time
Ts = 1;

% NOTE: this is a previous version of the system with lower
% coupling (-0.2) and epsilon = [0.01; 0.01].

% Discrete time state space model
A = [
      0.8890       0     1 -0.2
           0  0.8890  -0.2    1
           0       0     1    0
           0       0     0    1
];
B = [
         1 -0.2  0  0;
      -0.2    1  0  0;
         0    0  1  0;
         0    0  0  1
];
C = [ 0.1110 0         0  0;
             0  0.1110 0  0];
D = zeros(2, 4);
Gpss = ss(A,B,C,D,Ts);

% Dimensions
[n, nu, ny] = check_dimensions(A, B, C, D);

% Simulation settings
nT = 50;
t = Ts*(0:nT)';

% Choose time and amplitude of input disturbance
t_shock = [5 15];
du0 = [1; 1];
% When you make the shock larger the MKF observers
% do better
%du0 = [2; 2];

% Measured input
%U = (idinput(size(t)) + 1)/2;
U = zeros(nT+1, 2);
U(t >= 1, 1) = -1;

% Disturbance input
% This is used by the SKF observer
alpha = zeros(nT+1, 2);
alpha(t == t_shock(1), 1) = 1;
alpha(t == t_shock(2), 2) = 1;
Wp = du0' .* alpha;

U_sim = [U Wp];
x0 = zeros(n, 1);

% System models
sys_model.A = A;
sys_model.B = B;  % input disturbances unmeasured
sys_model.C = C;
sys_model.Ts = Ts;
sys_models = repmat({sys_model}, 1, 3);

% Choose measurement noise for plant
sigma_MP = [0; 0];  % See below for simulation with noise
V = sigma_MP'.*randn(nT+1, ny);

% Run simulation
[X, Y, Ym] = run_simulation_sys(sys_models,U_sim,V,alpha,nT,x0);

% % Simulate system
% X2 = zeros(nT+1,n);
% Y2 = zeros(nT+1,ny);
% xk = x0;
% for i = 1:nT+1
% 
%     % Inputs
%     uk = U_sim(i,:)';
% 
%     % Compute y(k)
%     yk = C * xk + D * uk;
% 
%     % Store results
%     X2(i, :) = xk';
%     Y2(i, :) = yk';
% 
%     % Compute x(k+1)
%     xk = A * xk + B * uk;
% 
% end
% 
% % Check simulation output is correct
% [Y3, t, X3] = lsim(Gpss, U_sim, t, x0);
% assert(isequal(X, X2))
% assert(isequal(Y, Y2))
% assert(isequal(X, X3))
% assert(isequal(Y, Y3))

% % Plot of inputs and outputs
% figure(4); clf
% 
% ax1 = subplot(5,1,1:2);
% plot(t,Y,'Linewidth',2); hold on
% plot(t,Ym,'.');
% max_min = [min(min([Y Ym])) max(max([Y Ym]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% ylabel('$y_i(k)$', 'Interpreter', 'latex')
% legend(compose("$y_%d(k)$", 1:ny), 'Interpreter', 'latex')
% title('System output and output measurements')
% grid on
% 
% P = cumsum(Wp);
% 
% ax2 = subplot(5,1,3:4);
% stairs(t,U,'Linewidth',2); hold on
% stairs(t,P,'Linewidth',2);
% max_min = [min(min([U P])) max(max([U P]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% ylabel('$u_i(k)$ and $w_{p_i}(k)$', 'Interpreter', 'latex')
% labels = [compose("$u_%d(k)$", 1:size(U,2)) ...
%     compose("$w_{p,%d}(k)$", 1:size(alpha,2))];
% legend(labels, 'Interpreter', 'latex')
% title('Input')
% grid on
% 
% ax3 = subplot(5,1,5);
% stairs(t,alpha,'Linewidth',2); hold on
% max_min = [min(min(alpha)) max(max(alpha))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$\gamma_i(k)$')
% legend(compose("$%s_%d(k)$", "\gamma", 1:size(alpha,2)), 'Interpreter', 'latex')
% title('Model sequence')
% grid on
% 
% linkaxes([ax1 ax2 ax3], 'x')

% Designate measured input and output signals
u_meas = [true; true; false; false];
y_meas = [true; true];

% Observer model (without unmeasured disturbance input)
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);
nw = sum(~u_meas);

% Disturbance input (used by SKF observer)
Bw = B(:, ~u_meas);
nw = sum(~u_meas);

% Dimensions of observer model
[n, nu, ny] = check_dimensions(A, Bu, C, Du);

% RODD random variable parameters
epsilon = [0.01; 0.01];
sigma_M = [0.05; 0.05];
%sigma_M = [0; 0];  % set to zero for testing
sigma_wp = [0.01 1; 0.01 1];

% Observer models
model.A = A;
model.B = Bu;  % input disturbances unmeasured
model.C = C;
model.Ts = Ts;
models = repmat({model}, 1, 3);

% Observer parameters (same for all observers)
models{1}.Q = diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,1)^2]);
models{2}.Q = diag([0.01 0.01 sigma_wp(1,2)^2 sigma_wp(2,1)^2]);
models{3}.Q = diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,2)^2]);
assert(isequal(size(models{1}.Q), [n n]))
assert(isequal(size(models{2}.Q), [n n]))
assert(isequal(size(models{3}.Q), [n n]))

R = diag(sigma_M.^2);
models{1}.R = R;
models{2}.R = R;
models{3}.R = R;
assert(isequal(size(R), [ny ny]))

P0 = 1000*eye(n);
%P0_init = repmat({P0}, 1, 3);
x0 = zeros(n,1);
y0 = models{1}.C * x0;
r0 = [1 1]';  % initial system mode

% Custom MKF test observers
% Devise a custom multi-model filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% this test simulation (t = t_shock)

% Multiple model filter 1
seq = repmat({ones(1, nT+1)}, 4, 1);
seq{2}(t == t_shock(1)) = 2;  % shock 1
seq{3}(t == t_shock(2)) = 3;  % shock 2
seq{4}(t == t_shock(1)) = 2;  % both
seq{4}(t == t_shock(2)) = 3;

% seq{2}(t == t_shock(1)) = 2;
% seq{2}(t == t_shock(2)) = 3;

% Build probability transition matrix
p_rk_g_rkm1 = [1-epsilon epsilon]';
Z = [1 1; 2 1; 1 2];  % combinations
p_rk_g_rkm1 = prod(prob_rk(Z', p_rk_g_rkm1), 1)';
p_rk_g_rkm1 = p_rk_g_rkm1 ./ sum(p_rk_g_rkm1);  % normalized
T = repmat(p_rk_g_rkm1', 3, 1);

MKF3 = MKFObserverS(models,P0,seq,T,'MKF3');
assert(MKF3.nh == size(seq, 1))

% Multiple model filter 2
seq = {ones(1, nT+1)};
seq{1}(t == t_shock(1)) = 2;
seq{1}(t == t_shock(2)) = 3;
MKF4 = MKFObserverS(models,P0,seq,T,'MKF4');
assert(MKF4.nh == 1)

% Define scheduled Kalman filter
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
% seq = sum(alpha .* 2.^(1:-1:0), 2)';
SKF = SKFObserverS(models,P0,seq{1},"SKF");

% Choose observers to test
observers = {MKF3, MKF4, SKF};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Identify which observer to log MKF data for
f_mkf = 1;

% Simulate observers
[Xk_est,Yk_est,DiagP,MKF_vars] = ...
    run_simulation_obs(Ym,U,alpha,seq,observers,f_mkf);

% Output estimation errors
E_obs = repmat(Y, 1, n_obs) - Yk_est;

% Combine and display results
%sim_results = table(t,alpha,U,X,Y,Ym,Xk_est,Yk_est,E_obs);
%disp(sim_results)

% % Plot observer estimates
% figure(5); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,obs_labels)

% % Plot MKF observer variables
% figure(6); clf
% 
% ax1 = subplot(2,1,1);
% ax_labels = {'Time ($t$)', 'Hypotheses, $f$', ...
%     '$\mathrm{Tr}(\mathbf{P}(k))$'};
% make_waterfall_plot(t, MKF_vars.trP_obs, [0 5], ax_labels, [0 64]);
% title("Trace of error covariance")
% 
% ax2 = subplot(2,1,2);
% ax_labels = {'Time ($t$)', 'Hypotheses, $f$', ...
%     '$\mathrm{Pr}(\Gamma_f(k) | \mathbf{Y}(k))$'};
% make_waterfall_plot(t, MKF_vars.p_seq_g_Yk, [0 1], ax_labels, [0 64]);
% title("Hypothesis probabilities")
% 
% %linkaxes([ax1 ax2], 'x')

MKF_Yk_est = cell2mat(cellfun(@(x) x', MKF_vars.Yk_est, ...
    'UniformOutput', false));

% % Plot MKF observer variables
% figure(7); clf
% nh = observers{f_mkf}.nh;
% 
% ax1 = subplot(2,1,1);
% plot(t, Y(:, 1), '--'); hold on
% plot(t, MKF_Yk_est(:, 1:2:nh*2));
% plot(t, Yk_est(:, 1), 'k.-');
% grid on
% ylabel("$y_1(k)$", 'Interpreter', 'latex')
% legend(["$y_1(k)$" compose("$%s{y}_{%d,1}(k)$", "\hat", 1:nh) ...
%     sprintf("$%s{y_1}(k)$", "\hat")], 'Interpreter', 'latex')
% 
% ax2 = subplot(2,1,2);
% plot(t, Y(:, 2), '--'); hold on
% plot(t, MKF_Yk_est(:, (1:2:nh*2) + 1));
% plot(t, Yk_est(:, 2), 'k.-');
% grid on
% ylabel("$y_2(k)$", 'Interpreter', 'latex')
% legend(["$y_2(k)$" compose("$%s{y}_{%d,2}(k)$", "\hat", 1:nh) ...
%     sprintf("$%s{y_2}(k)$", "\hat")], 'Interpreter', 'latex')
% xlabel('Time ($t$)', 'Interpreter', 'latex')
% 
% linkaxes([ax1 ax2], 'x')

% Compute mean-squared errors
MSE = struct();
for i = 1:n_obs
    MSE.(observers{i}.label) = mean(E_obs(:, i*ny-1:i*ny).^2);
end
%disp(MSE)

% Check final state estimates of all observers
test_X_est = [
   -1.7921    8.8475    0.9998    0.9997 ...
   -1.7921    8.8475    0.9998    0.9997 ...
   -1.7921    8.8475    0.9998    0.9997 ...
];
assert(isequal(round(Xk_est(t == t(end), :), 4), test_X_est))

% Check final error covariance estimates
test_DiagP = [ ...
    0.0405    0.0407    0.0015    0.0015 ...
    0.0405    0.0407    0.0015    0.0015 ...
    0.0405    0.0407    0.0015    0.0015 ...
];
assert(isequal(round(DiagP(t == t(end), :), 4), test_DiagP))

% Display trace of covariance matrix data for MKF observer filters
% disp([table(t) array2table(MKF_vars.trP_obs, 'VariableNames', ...
%     compose("Tr(P_%d)", 1:observers{f_mkf}.nh))])

% Results on 2022-09-08
MSE_test_values = struct( ...
 'MKF3', [0.000228, 0.000194],  ...
 'MKF4', [0.000009, 0.000010],  ...
 'SKF', [0.000009, 0.000010]  ...
);

labels = fieldnames(MSE);
for i = 1:numel(labels)
    fprintf("%s: %.6f, %.6f (%.6f, %.6f)\n", labels{i}, MSE.(labels{i}), ...
        MSE_test_values.(labels{i}))
end
for i = 1:numel(labels)
    assert(isequal(round(MSE.(labels{i}), 6), MSE_test_values.(labels{i})))
end

% Now re-run simulation with noise

% Choose measurement noise for plant
sigma_MP = [0.05; 0.05];  % See above for simulation without noise
V = sigma_MP'.*randn(nT+1, ny);

% Run simulation
[X, Y, Ym] = run_simulation_sys(sys_models,U_sim,V,alpha,nT,x0);

% Reset observers to initial states
for f = 1:n_obs
    observers{f}.reset()
end

% Simulate observers
[Xk_est,Yk_est,DiagP,MKF_vars] = ...
    run_simulation_obs(Ym,U,alpha,seq,observers,f_mkf);

% % Plot of inputs and outputs
% figure(8); clf
% 
% ax1 = subplot(5,1,1:2);
% plot(t,Y,'Linewidth',2); hold on
% plot(t,Ym,'.');
% max_min = [min(min([Y Ym])) max(max([Y Ym]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% ylabel('$y_i(k)$', 'Interpreter', 'latex')
% legend(compose("$y_%d(k)$", 1:ny), 'Interpreter', 'latex')
% title('System output and output measurements')
% grid on
% 
% P = cumsum(Wp);
% 
% ax2 = subplot(5,1,3:4);
% stairs(t,U,'Linewidth',2); hold on
% stairs(t,P,'Linewidth',2);
% max_min = [min(min([U P])) max(max([U P]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% ylabel('$u_i(k)$ and $w_{p_i}(k)$', 'Interpreter', 'latex')
% labels = [compose("$u_%d(k)$", 1:size(U,2)) ...
%     compose("$w_{p,%d}(k)$", 1:size(alpha,2))];
% legend(labels, 'Interpreter', 'latex')
% title('Input')
% grid on
% 
% ax3 = subplot(5,1,5);
% stairs(t,alpha,'Linewidth',2); hold on
% max_min = [min(min(alpha)) max(max(alpha))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$\gamma_i(k)$')
% legend(compose("$%s_%d(k)$", "\gamma", 1:size(alpha,2)), 'Interpreter', 'latex')
% title('Model sequence')
% grid on
% 
% linkaxes([ax1 ax2 ax3], 'x')

% Output estimation errors
E_obs = repmat(Y, 1, n_obs) - Yk_est;

% Combine and display results
sim_results = table(t,alpha,U,X,Y,Ym,Xk_est,Yk_est,E_obs);
%disp(sim_results)

% % Plot observer estimates
% figure(9); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,obs_labels)

% % Plot MKF observer variables
% figure(10); clf
% 
% ax1 = subplot(2,1,1);
% ax_labels = {'Time ($t$)', 'Hypotheses, $f$', ...
%     '$\mathrm{Tr}(\mathbf{P}(k))$'};
% make_waterfall_plot(t, MKF_vars.trP_obs, [0 5], ax_labels, [0 64]);
% title("Trace of error covariance")
% 
% ax2 = subplot(2,1,2);
% ax_labels = {'Time ($t$)', 'Hypotheses, $f$', ...
%     '$\mathrm{Pr}(\Gamma_f(k) | \mathbf{Y}(k))$'};
% make_waterfall_plot(t, MKF_vars.p_seq_g_Yk, [0 1], ax_labels, [0 64]);
% title("Hypothesis probabilities")
% 
% %linkaxes([ax1 ax2], 'x')

MKF_Yk_est = cell2mat(cellfun(@(x) x', MKF_vars.Yk_est, ...
    'UniformOutput', false));

% % Plot MKF observer variables
% figure(11); clf
% nh = observers{f_mkf}.nh;
% 
% ax1 = subplot(2,1,1);
% plot(t, Y(:, 1), '--'); hold on
% plot(t, MKF_Yk_est(:, 1:2:nh*2));
% plot(t, Yk_est(:, 1), 'k.-');
% grid on
% ylabel("$y_1(k)$", 'Interpreter', 'latex')
% legend(["$y_1(k)$" compose("$%s{y}_{%d,1}(k)$", "\hat", 1:nh) ...
%     sprintf("$%s{y_1}(k)$", "\hat")], 'Interpreter', 'latex')
% 
% ax2 = subplot(2,1,2);
% plot(t, Y(:, 2), '--'); hold on
% plot(t, MKF_Yk_est(:, (1:2:nh*2) + 1));
% plot(t, Yk_est(:, 2), 'k.-');
% grid on
% ylabel("$y_2(k)$", 'Interpreter', 'latex')
% legend(["$y_2(k)$" compose("$%s{y}_{%d,2}(k)$", "\hat", 1:nh) ...
%     sprintf("$%s{y_2}(k)$", "\hat")], 'Interpreter', 'latex')
% xlabel('Time ($t$)', 'Interpreter', 'latex')
% 
% linkaxes([ax1 ax2], 'x')

% Compute mean-squared errors
MSE = struct();
for i = 1:n_obs
    MSE.(observers{i}.label) = mean(E_obs(:, i*ny-1:i*ny).^2);
end
disp(MSE);

% Check final state estimates are close to true system values
assert(all(abs(repmat(X(end, :), 1, n_obs) - Xk_est(end, :)) < 0.2))

% Check final error covariance estimates
test_DiagP = [ ...
    0.0405    0.0407    0.0015    0.0015 ...
    0.0405    0.0407    0.0015    0.0015 ...
    0.0405    0.0407    0.0015    0.0015 ...
];
assert(isequal(round(DiagP(t == t(end), :), 4), test_DiagP))

% Display trace of covariance matrix data for MKF observer filters
% disp([table(t) array2table(MKF_vars.trP_obs, 'VariableNames', ...
%     compose("Tr(P_%d)", 1:observers{f_mkf}.nh))])

% Results on 2022-09-08
% Note: these results may depend somewhat on the random initialization
MSE_test_values = struct( ...
 'MKF3', [0.001051, 0.000828],  ...
 'MKF4', [0.000606, 0.000692],  ...
 'SKF', [0.000606, 0.000692]  ...
);

labels = fieldnames(MSE);
% for i = 1:numel(labels)
%     fprintf("%s: %.6f, %.6f (%.6f, %.6f)\n", labels{i}, MSE.(labels{i}), ...
%         MSE_test_values.(labels{i}))
% end
for i = 1:numel(labels)
    assert(isequal(round(MSE.(labels{i}), 6), MSE_test_values.(labels{i})))
end

