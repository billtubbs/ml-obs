% Tests the following observer classes
% 
%  - KalmanFilterF
%      Single Kalman Filter (used for comparison).
%  - SKFObserver
%      KF with switching system model.
%  - SKFObserverS
%      KF with switching system model and pre-set sequence.
%  - MKFObserver
%      Multiple-model KF with switching system models.
%  - MKFObserverS
%      Multiple-model KF with switching system models and
%      pre-set sequence.
%


clear all

seed = 0;
rng(seed)


%% Simulation test - SISO system

% Load switching system
sys_js2_siso

% Check dimensions
assert(isequal(size(model1.A), [n n]))
assert(isequal(size(model1.B), [n nu]))
assert(isequal(size(model1.C), [ny n]))
assert(isequal(size(model2.A), [n n]))
assert(isequal(size(model2.B), [n nu]))
assert(isequal(size(model2.C), [ny n]))
assert(nj == 2)

% Input disturbance variance
%sigma_w = 0.1;
sigma_w = 0;

% Process noise std. dev.
sigma_W = [0; 0];

% Measurement noise std. dev.
sigma_M = 0.0;

% Simulation settings
nT = 60;
t = Ts*(0:nT)';

% Inputs
%U = (idinput(size(t)) + 1)/2;
U = zeros(nT+1,1);
U(t>2) = 1;
V = sigma_M * randn(size(t));

% Switching sequence
%Gamma = int8(rand(nT+1, 1) > T(1, 1));
Gamma = int8(zeros(nT+1, 1));
Gamma(t>=10, 1) = 1;

% Simulate switching system
[X, Y, Ym] = run_simulation_sys(models,U,V,Gamma,nT);

% % Plot of inputs and outputs
% figure(1); clf
% 
% ax1 = subplot(5,1,1:2);
% plot(t,Y,'Linewidth',2); hold on
% plot(t,Ym,'o');
% max_min = [min(min([Y Ym])) max(max([Y Ym]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% ylabel('y(k)')
% title('System output and output measurements')
% grid on
% 
% ax2 = subplot(5,1,3:4);
% stairs(t,U,'Linewidth',2);
% max_min = [min(min(U)) max(max(U))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% ylabel('u(k) and w_p(k)')
% legend('u(k)')
% title('Input')
% grid on
% 
% ax3 = subplot(5,1,5);
% stairs(t,Gamma,'Linewidth',2)
% max_min = [min(min(Gamma)) max(max(Gamma))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('gamma(k)')
% title('Model sequence')
% grid on
% 
% linkaxes([ax1 ax2 ax3], 'x')

% Observer parameters (same for all observers)
P0 = 10000;
x0 = 0.5;
y0 = models{1}.C * x0;
r0 = 1;  % system mode
models{1}.Q = 0.01;
models{1}.R = 0.1^2;
models{2}.Q = 0.01;
models{2}.R = 0.1^2;
assert(isequal(size(models{1}.Q), size(models{2}.Q)))
assert(isequal(size(models{1}.R), size(models{2}.R)))

% Standard Kalman filters
KF1 = KalmanFilterF(models{1},P0,'KF1');
KF2 = KalmanFilterF(models{2},P0,'KF2');

% Transition probabilities
epsilon = 0.05;
T = [1-epsilon epsilon; epsilon 1-epsilon];
assert(all(sum(T, 2) == 1))

% System indicator sequences
seq1 = {
    ones(1, nT+1);
    [ones(1, 20) 2*ones(1, nT+1-20)];  % matches Gamma'
    [ones(1, 40) 2*ones(1, nT+1-40)];
    2*ones(1, nT+1);
 };
assert(isequal(seq1{2}, Gamma' + 1))

% Define observers with a switching system

% Define switching Kalman filter
SKF = SKFObserver(models,P0,"SKF1");

% Test initialisation
assert(strcmp(SKF.type, "SKF"))
assert(isequal(SKF.models, models))
assert(isequal(SKF.P0, P0))
assert(strcmp(SKF.label, "SKF1"))
assert(isequal(SKF.x0, zeros(n, 1)))
assert(isequal(SKF.r0, r0))
assert(SKF.nu == nu)
assert(SKF.ny == ny)
assert(SKF.nj == nj)
assert(SKF.n == n)
assert(isequal(SKF.xkp1_est, zeros(n, 1)))
assert(isequal(SKF.Pkp1, P0))
assert(isequal(SKF.rk, r0))
assert(isequal(SKF.rk, r0))
assert(isequaln(SKF.xk_est, nan(n, 1)))
assert(isequaln(SKF.Pk, nan(n)))
assert(isequaln(SKF.yk_est, nan(ny, 1)))
assert(isequaln(SKF.Kf, nan(n, ny)))
assert(isequaln(SKF.Sk, nan(ny)))

% Redefine this time with initial conditions
SKF = SKFObserver(models,P0,"SKF1",x0);
assert(isequal(SKF.xkp1_est, x0))

SKF = SKFObserver(models,P0,"SKF1",x0,r0);
assert(isequal(SKF.x0, x0))
assert(isequal(SKF.xkp1_est, x0))
assert(isequal(SKF.r0, r0))
assert(isequal(SKF.rk, r0))

% Again with default initial conditions
SKF = SKFObserver(models,P0,"SKF");

% Define scheduled SKF filter
seq = Gamma' + 1;  % add one for MATLAB indexing
SKF_S = SKFObserverS(models,P0,seq,"SKF_S");

assert(strcmp(SKF_S.type, "SKF_S"))
assert(isequal(SKF_S.models, models))
assert(isequal(SKF_S.P0, P0))
assert(strcmp(SKF_S.label, "SKF_S"))
assert(SKF_S.n == n)
assert(SKF_S.nu == nu)
assert(SKF_S.ny == ny)
assert(SKF_S.nj == nj)
assert(isequal(SKF_S.xkp1_est, zeros(n, 1)))
assert(isequal(SKF_S.Pkp1, P0))
assert(isequal(SKF_S.rk, seq(:, 1)))
assert(isequaln(SKF_S.xk_est, nan(n, 1)))
assert(isequaln(SKF_S.Pk, nan(n)))
assert(isequaln(SKF_S.yk_est, nan(ny, 1)))
assert(isequaln(SKF_S.Kf, nan(n, ny)))
assert(isequaln(SKF_S.Sk, nan(ny)))
assert(isequaln(SKF_S.seq, seq))
assert(isequaln(SKF_S.nf, size(seq, 2)))
assert(isequaln(SKF_S.i, 0))
assert(isequaln(SKF_S.i_next, 1))

% Transition probabilities
epsilon = 0.05;
T = [1-epsilon epsilon; epsilon 1-epsilon];
assert(all(sum(T, 2) == 1))

% System indicator sequences
seq1 = {
    ones(1, nT+1);
    [ones(1, 20) 2*ones(1, nT+1-20)];  % matches Gamma'
    [ones(1, 40) 2*ones(1, nT+1-40)];
    2*ones(1, nT+1);
 };
assert(isequal(seq1{2}, Gamma' + 1))

% Define MKF observer 1
seq = seq1;
nh = numel(seq);
%P0j = repmat({P0}, n_filt, 1);

% First, define with no initial state specified (should be set to zero)
% TODO: Allow independent P0 to be specified for each filter.
r0 = cellfun(@(s) s(:, 1), seq1);
MKF = MKFObserver(models,P0,T,r0,"MKF1");

% Test initialisation
assert(strcmp(MKF.type, "MKF"))
assert(isequal(MKF.models, models))
assert(isequal(MKF.P0, P0))
assert(isequal(MKF.Pkp1, P0))
assert(strcmp(MKF.label, "MKF1"))
assert(MKF.n == n)
assert(MKF.nu == nu)
assert(MKF.ny == ny)
assert(MKF.nj == nj)
assert(isequal(MKF.xkp1_est, zeros(n, 1)))
assert(MKF.nh == nh)
assert(isequal(MKF.p_seq_g_Yk_init, ones(nh, 1) ./ nh))
assert(isequal(MKF.rk, r0))
assert(isequaln(MKF.p_yk_g_seq_Ykm1, nan(nh, 1)))
assert(isequaln(MKF.p_rk_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF.p_seq_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF.filters.Xkp1_est, zeros(n, 1, nh)))
assert(isequaln(MKF.filters.Pkp1, repmat(P0, 1, 1, nh)))
assert(isequaln(MKF.filters.Kf, nan(n, ny, nh)))
assert(isequaln(MKF.filters.Sk, repmat(nan(ny), 1, 1, nh)))
assert(isequaln(MKF.xk_est, nan(n, 1)))
assert(isequaln(MKF.Pk, nan(n)))
assert(isequaln(MKF.yk_est, nan(ny, 1)))

% Redefine this time with initial conditions
MKF = MKFObserver(models,P0,T,r0,"MKF1",x0);
assert(isequal(MKF.x0, x0))
assert(isequal(MKF.r0, r0))
assert(isequal(MKF.xkp1_est, x0))
assert(isequal(MKF.rk, r0))
assert(isequaln(MKF.filters.Xkp1_est, repmat(x0, 1, 1, nh)))
assert(isequal(MKF.p_seq_g_Yk_init, ones(nh, 1) ./ nh))
assert(isequal(MKF.p_seq_g_Yk, MKF.p_seq_g_Yk_init))

% Test with initial hypothesis probabilities
p_seq_g_Yk_init = [0.6; 0.4];
MKF = MKFObserver(models,P0,T,r0,"MKF1",x0,p_seq_g_Yk_init);
assert(isequal(MKF.x0, x0))
assert(isequal(MKF.r0, r0))
assert(isequal(MKF.xkp1_est, x0))
assert(isequal(MKF.rk, r0))
assert(isequal(MKF.p_seq_g_Yk_init, p_seq_g_Yk_init))
assert(isequal(MKF.p_seq_g_Yk, p_seq_g_Yk_init))

% Again with default initial conditions
MKF = MKFObserver(models,P0,T,r0,"MKF");

% First, define with no initial state specified (should be set to zero)
% TODO: Allow independent P0 to be specified for each filter.
r0 = cellfun(@(s) s(:, 1), seq1);
MKF_S = MKFObserverS(models,P0,seq,T,"MKF_S");

% Test initialisation
assert(strcmp(MKF_S.type, "MKF-S"))
assert(isequal(MKF_S.models, models))
assert(isequal(MKF_S.P0, P0))
assert(isequal(MKF_S.Pkp1, P0))
assert(strcmp(MKF_S.label, "MKF_S"))
assert(MKF_S.n == n)
assert(MKF_S.nu == nu)
assert(MKF_S.ny == ny)
assert(MKF_S.nj == nj)
assert(isequal(MKF_S.xkp1_est, zeros(n, 1)))
assert(MKF_S.nh == nh)
assert(isequal(MKF_S.p_seq_g_Yk_init, ones(nh, 1) ./ nh))
assert(isequal(MKF_S.rk, r0))
assert(isequaln(MKF_S.p_yk_g_seq_Ykm1, nan(nh, 1)))
assert(isequaln(MKF_S.p_rk_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF_S.p_seq_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF_S.filters.Xkp1_est, zeros(n, 1, nh)))
assert(isequaln(MKF_S.filters.Pkp1, repmat(P0, 1, 1, nh)))
assert(isequaln(MKF_S.filters.Kf, nan(n, ny, nh)))
assert(isequaln(MKF_S.filters.Sk, repmat(nan(ny), 1, 1, nh)))
assert(isequaln(MKF_S.xk_est, nan(n, 1)))
assert(isequaln(MKF_S.Pk, nan(n)))
assert(isequaln(MKF_S.yk_est, nan(ny, 1)))
assert(isequaln(MKF_S.seq, seq))
assert(isequaln(MKF_S.nf, size(seq{1}, 2)))
assert(isequaln(MKF_S.i, 0))
assert(isequaln(MKF_S.i_next, 1))

% Define autonomous multi-model (AMM) observer 
MKF_AMM = MKFObserverAMM(models,P0,T,"MKF_AMM");

% Test initialisation
assert(strcmp(MKF_AMM.type, "MKF_AMM"))
assert(isequal(MKF_AMM.models, models))
assert(isequal(MKF_AMM.P0, P0))
assert(isequal(MKF_AMM.Pkp1, P0))
assert(strcmp(MKF_AMM.label, "MKF_AMM"))
assert(MKF_AMM.n == n)
assert(MKF_AMM.nu == nu)
assert(MKF_AMM.ny == ny)
assert(MKF_AMM.nj == nj)
assert(isequal(MKF_AMM.xkp1_est, zeros(n, 1)))
assert(MKF_AMM.nh == nj)
assert(isequal(MKF_AMM.p_seq_g_Yk_init, ones(nj, 1) ./ nj))
assert(isequal(MKF_AMM.rk, [1 2]'))
assert(isequaln(MKF_AMM.p_yk_g_seq_Ykm1, nan(nj, 1)))
assert(isequaln(MKF_AMM.p_rk_g_Ykm1, nan(nj, 1)))
assert(isequaln(MKF_AMM.p_seq_g_Ykm1, nan(nj, 1)))
assert(isequaln(MKF_AMM.filters.Xkp1_est, zeros(n, 1, nj)))
assert(isequaln(MKF_AMM.filters.Pkp1, repmat(P0, 1, 1, nj)))
assert(isequaln(MKF_AMM.filters.Kf, nan(n, ny, nj)))
assert(isequaln(MKF_AMM.filters.Sk, repmat(nan(ny), 1, 1, nj)))
assert(isequaln(MKF_AMM.xk_est, nan(n, 1)))
assert(isequaln(MKF_AMM.Pk, nan(n)))
assert(isequaln(MKF_AMM.yk_est, nan(ny, 1)))

% Choose observers to include in simulation
observers = {KF1, KF2, SKF, SKF_S, MKF, MKF_S, MKF_AMM};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Identify which observer to log MKF data for
f_mkf = [];

% Simulate observers - without measurement noise (Y)
[Xk_est,Yk_est,DiagP,MKF_K_obs,MKF_trP_obs,MKF_i,MKF_p_seq_g_Yk] = ...
    run_simulation_obs(Ym,U,Gamma,seq,observers,f_mkf);

% Output estimation errors
E_obs = Y - Yk_est;

% Combine and display results
sim_results1 = table(t,Gamma,U,X,Y,Ym,Xk_est,Yk_est,E_obs);

% figure(2); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,obs_labels)

% Check KF1 was accurate before system switched
assert(nanmean(E_obs(t < 10, 1).^2) < 0.0001)

% Check KF2 was accurate after system switched
assert(mean(E_obs(t > 12, 2).^2) < 0.001)

% Check which observers match KF1 before system switched
KF1_x_est = Xk_est(t == 9.5, 1);
assert(isequal(abs(Xk_est(t == 9.5, :) - KF1_x_est) < 0.0001, ...
    [true false true true true true true]))
KF1_diagP = sum(DiagP(t == 9.5, 1));
assert(isequal(abs(DiagP(t == 9.5, :) - KF1_diagP) < 0.0001, ...
    [true false true true true true true]))

% Check which observers match KF2 after system switched
KF2_x_est = Xk_est(t == 30, 2);
assert(isequal(abs(Xk_est(t == 30, :) - KF2_x_est) < 0.0001, ...
    [false true true true true true false]))
KF2_diagP = sum(DiagP(t == 30, 2));
assert(isequal(abs(DiagP(t == 30, :) - KF2_diagP) < 0.0001, ...
    [false true true true true true false]))

% Check MKF_AMM estimates match KF1 estimates
KF1_x_est = Xk_est(:, 1);
MKF_AMM_x_est = Xk_est(:, 7);
assert(all(abs(KF1_x_est - MKF_AMM_x_est) < 0.0001))

% Compute mean-squared error
mses = nanmean(E_obs.^2);
%disp(array2table(mses, 'VariableNames', obs_labels))

% Check MKF and SKF observer estimation errors
assert(isequal(round(mses, 6), [3.806151 0.269363 0 0 0 0 3.806151]))
% Previously, using ykp1_est (i.e. prior predictions),
% MSEs were: [5.1728 0.4313 0.1296 0.0660]

% Check final values
XYk_est = cell2mat(cellfun(@(obs) [obs.xk_est obs.yk_est], ...
    observers', 'UniformOutput', false));
assert(isequal(round(XYk_est, 6), [ ...
     -1.191082  -0.357325
      9.901224  -2.970367
      9.901227  -2.970368
      9.901227  -2.970368
      9.901227  -2.970368
      9.901227  -2.970368
     -1.191082  -0.357325 ...
]))

% Reset observer states to initial conditions
KF1.reset()
KF2.reset()
SKF.reset();
SKF_S.reset();
MKF.reset()
MKF_S.reset()
MKF_AMM.reset()

% Check main variables reset to initial conditions
XYk_est = cell2mat(cellfun(@(obs) [obs.xk_est' obs.yk_est'], ...
    observers', 'UniformOutput', false));
assert(isequaln(XYk_est, nan(n_obs, n+ny)))
Xkp1_est = cell2mat(cellfun(@(obs) obs.xkp1_est', observers', ...
    'UniformOutput', false));
assert(all(Xkp1_est == 0))
Pk = cell2mat(cellfun(@(obs) obs.Pk, observers', ...
    'UniformOutput', false));
assert(isequaln(Pk, nan(n_obs, 1)))
Pkp1 = cell2mat(cellfun(@(obs) obs.Pkp1, observers', ...
    'UniformOutput', false));
assert(all(Pkp1 == P0))

% Check other variables
assert(SKF.rk == 1)
assert(SKF_S.rk == 1)
assert(isequal(MKF.rk, [1 1 1 2]'))
assert(isequal(MKF_S.rk, [1 1 1 2]'))
assert(isequal(MKF_AMM.rk, [1 2]'))
assert(isequal([SKF_S.i SKF_S.i_next], [0 1]))
assert(isequal([MKF_S.i MKF_S.i_next], [0 1]))
assert(isequal(MKF.p_seq_g_Yk, [0.25 0.25 0.25 0.25]'))
assert(isequal(MKF_S.p_seq_g_Yk, [0.25 0.25 0.25 0.25]'))
assert(isequal(MKF_AMM.p_seq_g_Yk, [0.5 0.5]'))

% Test copy methods

% Redefine a new observer (identical to existing)
MKF_new = MKFObserver(models,P0,T,r0,"MKF");
assert(isequaln(MKF_new, MKF))
MKF_new.label = "MKF_new";

% Make a copy
MKF_copy = MKF_new.copy();
assert(isequaln(MKF_copy, MKF_new))
MKF_copy.label = "MKF_copy";

% Run simulation of original, new and copy

% Choose observers to include in simulation
observers = {MKF, MKF_new, MKF_copy};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

[Xk_est,Yk_est,DiagP,MKF_K_obs,MKF_trP_obs,MKF_i,MKF_p_seq_g_Yk] = ...
    run_simulation_obs(Ym,U,Gamma,seq,observers,f_mkf);

% Output estimation errors
E_obs = Y - Yk_est;

% Combine and display results
sim_results2 = table(t,Gamma,U,X,Y,Ym,Xk_est,Yk_est,E_obs);

% Check estimates are identical
assert(all(sim_results2.Yk_est(:, 1) == sim_results2.Yk_est, [1 2]))
assert(all(sim_results2.Xk_est(:, 1) == sim_results2.Xk_est, [1 2]))

% figure(3); clf
% plot_obs_estimates(t,X,X_est,Y,Y_est,obs_labels)


%% Simulation test on 2x2 system

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
nT = 200;
t = Ts*(0:nT)';

% Choose time and amplitude of input disturbance
t_shock = [5 10];
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
model.A = A;
model.B = B;  % input disturbances unmeasured
model.C = C;
model.Ts = Ts;
models = repmat({model}, 1, nj);

% Choose measurement noise for plant
sigma_MP = [0; 0];  % Set to zero for testing
V = sigma_MP'.*randn(nT+1, ny);

% Run simulation
[X, Y, Ym] = run_simulation_sys(models,U_sim,V,alpha,nT,x0);

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
sigma_M = [0.1; 0.1];
%sigma_M = [0; 0];  % set to zero for testing
sigma_wp = [0.01 1; 0.01 1];

% Observer models
model.A = A;
model.B = Bu;  % input disturbances unmeasured
model.C = C;
model.Ts = Ts;
nj = 3;
models = repmat({model}, 1, nj);

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
seq = {ones(1, nT+1); ones(1, nT+1); ones(1, nT+1); ones(1, nT+1)};
seq{2}(t == t_shock(1)) = 2;  % shock 1
seq{3}(t == t_shock(2)) = 3;  % shock 2
seq{4}(t == t_shock(1)) = 2;  % both
seq{4}(t == t_shock(2)) = 3;
p_rk = [1-epsilon epsilon]';
Z = [1 1; 2 1; 1 2];  % combinations
p_rk = prod(prob_rk(Z', p_rk), 1)';
p_rk = p_rk ./ sum(p_rk);  % normalized
T = repmat(p_rk', 3, 1);

MKF3 = MKFObserverS(models,P0,seq,T,'MKF3');
assert(MKF3.nh == 4)

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
[Xk_est,Yk_est,DiagP,MKF_K_obs,MKF_trP_obs,MKF_i,MKF_p_seq_g_Yk] = ...
    run_simulation_obs(Ym,U,alpha,seq,observers,f_mkf);

% Output estimation errors
E_obs = repmat(Y, 1, n_obs) - Yk_est;

% Combine and display results
sim_results = table(t,alpha,U,X,Y,Ym,Xk_est,Yk_est,E_obs);
disp(sim_results)

% Plot observer estimates
% figure(5); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,obs_labels)

% Check final state estimates
% Previous values (probably run with noise):
% test_X_est = [-1.801802  9.009008  1.000000  1.000000 -1.801802 ...
%     9.009008  1.000000  1.000000 -1.801802  9.009008  1.000000  1.000000];
test_X_est = [
    -1.801832  9.008979  0.999994  0.999994 ...
    -1.801802  9.009008  1.000000  1.000000 ...
    -1.801802  9.009008  1.000000  1.000000];
assert(isequal(round(Xk_est(t == t(end), :), 6), test_X_est))

% Check final error covariance estimates
% TODO: Haven't checked if these are correct.
% Previous values (probably run with noise):
% test_DiagP = [ 0.092947  0.092947  0.002086  0.002086  0.092947 ...
%     0.092947  0.002086  0.002086  0.092947  0.092947  0.002086  0.002086];
test_DiagP = [ ...
    0.083317    0.083317    0.001986    0.001986 ...
    0.083317    0.083317    0.001986    0.001986 ...
    0.083317    0.083317    0.001986    0.001986];
assert(isequal(round(DiagP(t == t(end), :), 6), test_DiagP))

% Display trace of covariance matrix data for MKF observer filters
[table(t) array2table(MKF_trP_obs, 'VariableNames', ...
    compose("Tr(P_%d)", 1:observers{f_mkf}.nh))]

% Compute mean-squared errors
MSE = struct();
for i = 1:n_obs
    MSE.(observers{i}.label) = mean(E_obs(:, i*ny-1:i*ny).^2);
end
disp(MSE)

% Results on Nov 8 after reverting back the Bayesian updating
MSE_test_values = struct( ...
 'MKF3', [0.000707 0.000348],  ...  % TODO: Why is it [0.001097 0.003716]
 'MKF4', [0.000017 0.000022],  ...
 'SKF', [0.000017 0.000022]  ...
);

labels = fieldnames(MSE);
for i = 1:numel(labels)
    fprintf("%s: %.6f, %.6f (%.6f, %.6f)\n", labels{i}, MSE.(labels{i}), ...
        MSE_test_values.(labels{i}))
end
for i = 1:numel(labels)
    assert(isequal(round(MSE.(labels{i}), 6), MSE_test_values.(labels{i})))
end


%% Test copy methods

% Load switching system
sys_js2_siso

% Observer parameters (same for all observers)
P0 = 10000;
x0 = 0.5;
r0 = 1;  % system mode
models{1}.Q = 0.01;
models{1}.R = 0.1^2;
models{2}.Q = 0.01;
models{2}.R = 0.1^2;

% Transition probabilities
epsilon = 0.05;
T = [1-epsilon epsilon; epsilon 1-epsilon];
assert(all(sum(T, 2) == 1))

% System indicator sequences
nT = 60;
seq1 = {
    ones(1, nT+1);
    [ones(1, 20) 2*ones(1, nT+1-20)];  % matches Gamma'
    [ones(1, 40) 2*ones(1, nT+1-40)];
    2*ones(1, nT+1);
 };

% Define MKF observer
seq = seq1;
n_filt = numel(seq);
%P0j = repmat({P0}, n_filt, 1);

% Define multi-model observer with initial conditions
MKF = MKFObserverS(models,P0,seq,T,'MKF',x0);

% Test handle copy
MKF_hcopy = MKF;
assert(isequaln(MKF_hcopy, MKF))  % same values
assert(MKF_hcopy == MKF)  % must be same object

MKF.x0 = 1.0;
assert(isequal(MKF_hcopy.x0, 1.0))

% Test true copy
MKF_copy = MKF.copy();
assert(isequaln(MKF_copy, MKF))  % same values
assert(MKF_copy ~= MKF)  % must not be same object
assert(isequaln(MKF_copy.filters, MKF.filters))
assert(isequaln(MKF_copy.filters, MKF.filters))

% Check deep copy was made
% TODO: This is not working
%assert(MKF_copy.filters{1} ~= MKF.filters{1})  % must not be same object
%assert(MKF_copy.filters{2} ~= MKF.filters{2})

MKF.label = "New name";
assert(~isequal(MKF_copy.label, "New name"))

MKF.filters.Xkp1_est(:, :, 1) = 99;
assert(~isequal(MKF_copy.filters.Xkp1_est(:, :, 1), 99))

%END


% TODO: Can't the function run_test_simulation.m be used here?
function [Xk_est,Yk_est,DiagP,MKF_K_obs,MKF_trP_obs,MKF_i,MKF_p_seq_g_Yk] = ...
    run_simulation_obs(Ym,U,Gamma,seq,observers,f_mkf)
% Simulate observers

    nT = size(Ym, 1) - 1;
    ny = size(Ym, 2);
    n_obs = numel(observers);
    n = size(observers{1}.xk_est, 1);

    Xk_est = zeros(nT+1, n*n_obs);
    Yk_est = zeros(nT+1, ny*n_obs);
    DiagP = zeros(nT+1, n*n_obs);

    if ~isempty(f_mkf)
        obs_mkf = observers{f_mkf};
        nh = size(obs_mkf.seq, 1);
        MKF_K_obs = cell(nT+1, n_obs*nh);
        MKF_trP_obs = nan(nT+1, nh);
        MKF_i = nan(nT+1, 2);
        MKF_p_seq_g_Yk = nan(nT+1, nh);
    else
        MKF_K_obs = {};
        MKF_trP_obs = nan;
        MKF_i = nan;
        MKF_p_seq_g_Yk = nan;
    end

    seq = cell2mat(seq);  % makes it easier to index

    for i = 1:nT+1

        yk = Ym(i, :)';
        uk = U(i, :)';

        % Update observers
        for f = 1:n_obs
            obs = observers{f};
            switch obs.type
                case "SKF"
                    rk = Gamma(i) + 1;
                    obs.update(yk, uk, rk);
                case "MKF"
                    rk = seq(:, i);
                    obs.update(yk, uk, rk);
                otherwise
                    obs.update(yk, uk);
            end
            if f == f_mkf
                MKF_i(i, :) = obs.i;
                MKF_p_seq_g_Yk(i, :) = obs.p_seq_g_Yk';
                for j = 1:obs.nh
                    MKF_K_obs{i, j} = obs.filters.Kf(:,:,j);
                    MKF_trP_obs(i, j) = trace(obs.filters.Pkp1(:,:,j));
                end
            end
            xk_est(1, (f-1)*n+1:f*n) = obs.xk_est';
            yk_est(1, (f-1)*ny+1:f*ny) = obs.yk_est';
            diagP(1, (f-1)*n+1:f*n) = diag(obs.Pk)';
        end

        % Record observer estimates
        Xk_est(i, :) = xk_est;
        Yk_est(i, :) = yk_est;
        DiagP(i, :) = diagP;

    end
end
