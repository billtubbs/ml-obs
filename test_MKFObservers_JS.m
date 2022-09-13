% Tests the following observer classes on a simple SISO jump system
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
%  - MKFObserverAMM
%      Autonomous (no switching) multiple-model KF with one
%      hypothesis for each system mode (i.e. nh = nj).
%  - MKFObserverGPB1
%      Suboptimal observer using generalised pseudo-Bayes 
%      algorithm - first order.
%


clear all

addpath("~/ml-plot-utils/")

seed = 0;
rng(seed)


%% Simulation test - SISO jump linear system

% Load system
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
assert(isempty(SKF.rkm1))
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
assert(isempty(SKF.rkm1))

% Finally, again with default initial conditions
SKF = SKFObserver(models,P0,"SKF");

% Define scheduled SKF filter
seq = Gamma' + 1;  % add one for MATLAB indexing
SKF_S = SKFObserverS(models,P0,seq,"SKF_S");

assert(strcmp(SKF_S.type, "SKF_S"))
assert(isequal(SKF_S.models, models))
assert(isequal(SKF_S.P0, P0))
assert(isequaln(SKF_S.seq, seq))
assert(strcmp(SKF_S.label, "SKF_S"))
assert(isequal(SKF.x0, 0))
assert(SKF_S.n == n)
assert(SKF_S.nu == nu)
assert(SKF_S.ny == ny)
assert(SKF_S.nj == nj)
assert(isequal(SKF_S.xkp1_est, zeros(n, 1)))
assert(isequal(SKF_S.Pkp1, P0))
assert(isequal(SKF_S.rk, seq(:, 1)))
assert(isempty(SKF_S.rkm1))
assert(isequaln(SKF_S.xk_est, nan(n, 1)))
assert(isequaln(SKF_S.Pk, nan(n)))
assert(isequaln(SKF_S.yk_est, nan(ny, 1)))
assert(isequaln(SKF_S.Kf, nan(n, ny)))
assert(isequaln(SKF_S.Sk, nan(ny)))
assert(isequaln(SKF_S.nf, size(seq, 2)))
assert(isequaln(SKF_S.i, 0))
assert(isequaln(SKF_S.i_next, 1))

% Transition probabilities
epsilon = 0.05;
T = [1-epsilon epsilon; epsilon 1-epsilon];
assert(all(sum(T, 2) == 1))

% System mode indicator sequences
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
assert(isempty(MKF.rkm1))
assert(isequaln(MKF.p_yk_g_seq_Ykm1, nan(nh, 1)))
assert(isequaln(MKF.p_rk_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF.p_rk_g_rkm1, nan(nh, 1)))
assert(isequaln(MKF.p_seq_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF.filters.Xkp1_est, zeros(n, 1, nh)))
assert(isequaln(MKF.filters.Pkp1, repmat(P0, 1, 1, nh)))
assert(isequaln(MKF.filters.Kf, nan(n, ny, nh)))
assert(isequaln(MKF.filters.Sk, nan(ny, ny, nh)))
assert(isequaln(MKF.xk_est, nan(n, 1)))
assert(isequaln(MKF.Pk, nan(n)))
assert(isequaln(MKF.yk_est, nan(ny, 1)))

% Redefine this time with initial conditions
MKF = MKFObserver(models,P0,T,r0,"MKF1",x0);
assert(isequal(MKF.x0, x0))
assert(isequal(MKF.r0, r0))
assert(isequal(MKF.xkp1_est, x0))
assert(isequal(MKF.rk, r0))
assert(isempty(MKF.rkm1))
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
assert(isempty(MKF.rkm1))
assert(isequal(MKF.p_seq_g_Yk_init, p_seq_g_Yk_init))
assert(isequal(MKF.p_seq_g_Yk, p_seq_g_Yk_init))

% Finally, again with default initial conditions
MKF = MKFObserver(models,P0,T,r0,"MKF");

% Define MKF observer 2 - with pre-determined switching sequence

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
assert(isequal(MKF_S.x0, zeros(n, 1)))
assert(MKF_S.n == n)
assert(MKF_S.nu == nu)
assert(MKF_S.ny == ny)
assert(MKF_S.nj == nj)
assert(isequal(MKF_S.xkp1_est, zeros(n, 1)))
assert(MKF_S.nh == nh)
assert(isequal(MKF_S.p_seq_g_Yk_init, ones(nh, 1) ./ nh))
assert(isequal(MKF_S.rk, r0))
assert(isempty(MKF_S.rkm1))
assert(isequaln(MKF_S.p_yk_g_seq_Ykm1, nan(nh, 1)))
assert(isequaln(MKF_S.p_rk_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF_S.p_rk_g_rkm1, nan(nh, 1)))
assert(isequaln(MKF_S.p_seq_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF_S.filters.Xkp1_est, zeros(n, 1, nh)))
assert(isequaln(MKF_S.filters.Pkp1, repmat(P0, 1, 1, nh)))
assert(isequaln(MKF_S.filters.Kf, nan(n, ny, nh)))
assert(isequaln(MKF_S.filters.Sk, nan(n, ny, nh)))
assert(isequaln(MKF_S.xk_est, nan(n, 1)))
assert(isequaln(MKF_S.Pk, nan(n)))
assert(isequaln(MKF_S.yk_est, nan(ny, 1)))
assert(isequaln(MKF_S.seq, seq))
assert(isequaln(MKF_S.nf, size(seq{1}, 2)))
assert(isequaln(MKF_S.i, 0))
assert(isequaln(MKF_S.i_next, 1))

% Define autonomous multi-model (AMM) observer 
MKF_AMM = MKFObserverAMM(models,P0,"MKF_AMM");

% Test initialisation
assert(strcmp(MKF_AMM.type, "MKF_AMM"))
assert(isequal(MKF_AMM.models, models))
assert(isequal(MKF_AMM.P0, P0))
assert(isequal(MKF_AMM.Pkp1, P0))
assert(strcmp(MKF_AMM.label, "MKF_AMM"))
assert(isequal(MKF_S.x0, zeros(n, 1)))
assert(MKF_AMM.n == n)
assert(MKF_AMM.nu == nu)
assert(MKF_AMM.ny == ny)
assert(MKF_AMM.nj == nj)
assert(isequal(MKF_AMM.xkp1_est, zeros(n, 1)))
assert(MKF_AMM.nh == nj)
assert(isequal(MKF_AMM.p_seq_g_Yk_init, ones(nj, 1) ./ nj))
assert(isequal(MKF_AMM.rk, [1 2]'))
assert(isempty(MKF_AMM.rkm1))
assert(isequaln(MKF_AMM.p_yk_g_seq_Ykm1, nan(nj, 1)))
assert(isequaln(MKF_AMM.p_rk_g_Ykm1, nan(nj, 1)))
assert(isequaln(MKF_AMM.p_rk_g_rkm1, nan(nj, 1)))
assert(isequaln(MKF_AMM.p_seq_g_Ykm1, nan(nj, 1)))
assert(isequaln(MKF_AMM.filters.Xkp1_est, zeros(n, 1, nj)))
assert(isequaln(MKF_AMM.filters.Pkp1, repmat(P0, 1, 1, nj)))
assert(isequaln(MKF_AMM.filters.Kf, nan(n, ny, nj)))
assert(isequaln(MKF_AMM.filters.Sk, nan(ny, 1, nj)))
assert(isequaln(MKF_AMM.xk_est, nan(n, 1)))
assert(isequaln(MKF_AMM.Pk, nan(n)))
assert(isequaln(MKF_AMM.yk_est, nan(ny, 1)))

% Define MKF observer with GPB1 algorithm

% First, define with no initial state specified (should be set to zero)
MKF_GPB1 = MKFObserverGPB1(models,P0,T,"MKF_GPB1");

% Test initialisation
assert(strcmp(MKF_GPB1.type, "MKF_GPB1"))
assert(isequal(MKF_GPB1.models, models))
assert(isequal(MKF_GPB1.P0, P0))
assert(isequal(MKF_GPB1.Pkp1, P0))
assert(strcmp(MKF_GPB1.label, "MKF_GPB1"))
assert(isequal(MKF_GPB1.x0, zeros(n, 1)))
assert(MKF_GPB1.n == n)
assert(MKF_GPB1.nu == nu)
assert(MKF_GPB1.ny == ny)
assert(MKF_GPB1.nj == nj)
assert(MKF_GPB1.nh == nj)
assert(isequal(MKF_GPB1.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_GPB1.p_seq_g_Yk_init, ones(nj, 1) ./ nj))
assert(isequal(MKF_GPB1.rk, [1 2]'))
assert(isempty(MKF_GPB1.rkm1))
assert(isequaln(MKF_GPB1.p_yk_g_seq_Ykm1, nan(nj, 1)))
assert(isequaln(MKF_GPB1.p_rk_g_Ykm1, nan(nj, 1)))
assert(isequaln(MKF_GPB1.p_rk_g_rkm1, nan(nj, 1)))
assert(isequaln(MKF_GPB1.p_seq_g_Ykm1, nan(nj, 1)))
assert(isequaln(MKF_GPB1.filters.Xkp1_est, zeros(n, 1, nj)))
assert(isequaln(MKF_GPB1.filters.Pkp1, repmat(P0, 1, 1, nj)))
assert(isequaln(MKF_GPB1.filters.Kf, nan(n, ny, nj)))
assert(isequaln(MKF_GPB1.filters.Sk, nan(ny, 1, nj)))
assert(isequaln(MKF_GPB1.xk_est, nan(n, 1)))
assert(isequaln(MKF_GPB1.Pk, nan(n)))
assert(isequaln(MKF_GPB1.yk_est, nan(ny, 1)))

% Define MKF observer with GPB2 algorithm

% First, define with no initial state specified (should be set to zero)
MKF_GPB2 = MKFObserverGPB2(models,P0,T,"MKF_GPB2");

% Test initialisation
assert(strcmp(MKF_GPB2.type, "MKF_GPB2"))
assert(isequal(MKF_GPB2.models, models))
assert(isequal(MKF_GPB2.P0, P0))
assert(isequal(MKF_GPB2.Pkp1, P0))
assert(strcmp(MKF_GPB2.label, "MKF_GPB2"))
assert(isequal(MKF_GPB2.x0, zeros(n, 1)))
assert(MKF_GPB2.n == n)
assert(MKF_GPB2.nu == nu)
assert(MKF_GPB2.ny == ny)
assert(MKF_GPB2.nj == nj)
assert(MKF_GPB2.nh == nj*nj)
assert(isequal(MKF_GPB2.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_GPB2.p_seq_g_Yk_init, ones(4, 1) ./ 4))
assert(isequal(MKF_GPB2.rk, [1 2 1 2]'))
assert(isequal(MKF_GPB2.rkm1, [1 1 2 2]'))
assert(isequaln(MKF_GPB2.p_yk_g_seq_Ykm1, nan(4, 1)))
assert(isequaln(MKF_GPB2.p_rk_g_Ykm1, nan(4, 1)))
assert(isequaln(MKF_GPB2.p_rk_g_rkm1, nan(4, 1)))
assert(isequaln(MKF_GPB2.p_seq_g_Ykm1, nan(4, 1)))
assert(isequaln(MKF_GPB2.filters.Xkp1_est, zeros(n, 1, 4)))
assert(isequaln(MKF_GPB2.filters.Pkp1, repmat(P0, 1, 1, 4)))
assert(isequaln(MKF_GPB2.filters.Kf, nan(n, ny, 4)))
assert(isequaln(MKF_GPB2.filters.Sk, nan(ny, 1, 4)))
assert(isequaln(MKF_GPB2.xk_est, nan(n, 1)))
assert(isequaln(MKF_GPB2.Pk, nan(n)))
assert(isequaln(MKF_GPB2.yk_est, nan(ny, 1)))

% Choose observers to include in simulation
observers = {KF1, KF2, SKF, SKF_S, MKF, MKF_S, MKF_AMM, MKF_GPB1, MKF_GPB2};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Identify which observer to log MKF data for
f_mkf = [];

% Simulate observers - without measurement noise (Y)
[Xk_est,Yk_est,DiagP,MKF_vars] = ...
    run_simulation_obs(Ym,U,Gamma,seq,observers,f_mkf);

% Output estimation errors
E_obs = Y - Yk_est;

% Combine and display results
sim_results1 = table(t,Gamma,U,X,Y,Ym,Xk_est,Yk_est,E_obs);

% Show in seperate plots so they are not too crowded
figure(2); clf
plot_obs_estimates(t,X,Xk_est(:, 1:4),Y,Yk_est(:, 1:4),obs_labels(1:4))
figure(3); clf
plot_obs_estimates(t,X,Xk_est(:, 5:end-1),Y,Yk_est(:, 5:end-1),obs_labels(5:end-1))
figure(4); clf
plot_obs_estimates(t,X,Xk_est(:, 7),Y,Yk_est(:, 7),obs_labels(7))
figure(5); clf
plot_obs_estimates(t,X,Xk_est(:, end),Y,Yk_est(:, end),obs_labels(end))

% Check KF1 was accurate before system switched
assert(nanmean(E_obs(t < 10, 1).^2) < 0.0001)

% Check KF2 was accurate after system switched
assert(mean(E_obs(t > 12, 2).^2) < 0.001)

% Check which observers match KF1 before system switched
KF1_x_est = Xk_est(t == 9.5, 1);
assert(isequal(abs(Xk_est(t == 9.5, :) - KF1_x_est) < 0.0001, ...
    [true false true true true true true true true]))
KF1_diagP = sum(DiagP(t == 9.5, 1));
assert(isequal(abs(DiagP(t == 9.5, :) - KF1_diagP) < 0.0001, ...
    [true false true true true true true true true]))

% Check which observers match KF2 after system switched
KF2_x_est = Xk_est(t == 30, 2);
assert(isequal(abs(Xk_est(t == 30, :) - KF2_x_est) < 0.0001, ...
    [false true true true true true false true true]))
KF2_diagP = sum(DiagP(t == 30, 2));
assert(isequal(abs(DiagP(t == 30, :) - KF2_diagP) < 0.0001, ...
    [false true true true true true false true true]))

% Check MKF_AMM estimates match KF1 estimates
KF1_x_est = Xk_est(:, 1);
MKF_AMM_x_est = Xk_est(:, 7);
assert(all(abs(KF1_x_est - MKF_AMM_x_est) < 0.0001))

% Compute mean-squared error
mses = nanmean(E_obs.^2);
%disp(array2table(mses, 'RowNames', {'MSE'}, 'VariableNames', obs_labels))

% Check MKF and SKF observer estimation errors
assert(isequal(round(mses, 6), ...
    [3.806151 0.269363 0 0 0 0 3.806151 0.000959 0.000959]))

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
     -1.191082  -0.357325
      9.901228  -2.970368
      9.901228  -2.970368 ...
]))

% Reset observer states to initial conditions
for f = 1:n_obs
    observers{f}.reset();
end

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
assert(isempty(SKF.rkm1))
assert(isempty(SKF_S.rkm1))
assert(isempty(MKF.rkm1))
assert(isempty(MKF_S.rkm1))
assert(isempty(MKF_AMM.rkm1))
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

% Run simulation of original, new, and copy

% Choose observers to include in simulation
observers = {MKF, MKF_new, MKF_copy};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

[Xk_est,Yk_est,DiagP,MKF_vars] = ...
    run_simulation_obs(Ym,U,Gamma,seq,observers,f_mkf);

% Output estimation errors
E_obs = Y - Yk_est;

% Combine and display results
sim_results2 = table(t,Gamma,U,X,Y,Ym,Xk_est,Yk_est,E_obs);

% Check estimates are identical
assert(all(sim_results2.Yk_est(:, 1) == sim_results2.Yk_est, [1 2]))
assert(all(sim_results2.Xk_est(:, 1) == sim_results2.Xk_est, [1 2]))

% figure(4); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,obs_labels)

return


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
% Not relevant when object does not contain handle objects
% assert(MKF_copy.filters ~= MKF.filters)  % must not be same object
% assert(MKF_copy.filters ~= MKF.filters)

MKF.label = "New name";
assert(~isequal(MKF_copy.label, "New name"))

MKF.filters.Xkp1_est(:, :, 1) = 99;
assert(~isequal(MKF_copy.filters.Xkp1_est(:, :, 1), 99))

