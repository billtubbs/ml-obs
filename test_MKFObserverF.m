% Test classes MKFObserverF and MKFObserverSchedF

clear all

seed = 0;
rng(seed)


%% Simulation test - SISO system

% Define system

% Sample period
Ts = 0.5;

% Discrete time state space models
% Model #1
A1 = 0.7;
B1 = 1;
C1 = 0.3;
D1 = 0;
Gpss1 = ss(A1,B1,C1,D1,Ts);

% Model #2
A2 = 0.9;
B2 = 1;
C2 = -0.3;  % -ve gain!
D2 = 0;
Gpss2 = ss(A2,B2,C2,D2,Ts);

% Dimensions
n = size(A1, 1);
nu = size(B1, 2);
ny = size(C1, 1);

% Check dimensions
assert(isequal(size(A1), size(A2)))
assert(isequal(size(B1), size(B2)))
assert(isequal(size(C1), size(C2)))
assert(isequal(size(D1), size(D2)))

% Define system models
m1.A = A1;
m1.B = B1;
m1.C = C1;
m1.Ts = Ts;
m2.A = A2;
m2.B = B2;
m2.C = C2;
m2.Ts = Ts;
models = {m1, m2};

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
% xlabel('t')
% ylabel('y(k)')
% title('System output and output measurements')
% grid on
% 
% ax2 = subplot(5,1,3:4);
% stairs(t,U,'Linewidth',2);
% max_min = [min(min(U)) max(max(U))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
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
models{1}.Q = 0.01;
models{1}.R = 0.1^2;
models{2}.Q = 0.01;
models{2}.R = 0.1^2;

% Standard Kalman filters
KF1 = KalmanFilterF(models{1},P0,'KF1',x0);
KF2 = KalmanFilterF(models{2},P0,'KF2',x0);

% Define scheduled MKF filter
seq = Gamma' + 1;
SKF1 = SKFObserverS(models,P0,seq,"SKF1",x0);
SKF2 = SKFObserverS(models,P0,seq,"SKF2",x0);

assert(strcmp(SKF1.type, "SKF_S"))
assert(isequal(SKF1.models, models))
assert(isequal(SKF1.Ts, Ts))
assert(isequal(SKF1.P0, P0))
assert(isequal(SKF1.Pkp1, P0))
assert(isequal(SKF1.seq, seq))
assert(strcmp(SKF1.label, "SKF1"))
assert(SKF1.n == n)
assert(SKF1.nu == nu)
assert(SKF1.ny == ny)
assert(SKF1.nf == size(SKF1.seq, 2))
assert(SKF1.nj == 2)
assert(isequal(SKF1.xkp1_est, x0))
assert(isequal(SKF1.rk, 1))

assert(strcmp(SKF2.type, "SKF_S"))
assert(isequal(SKF2.models, models))
assert(isequal(SKF2.Ts, Ts))
assert(isequal(SKF2.P0, P0))
assert(isequal(SKF2.seq, seq))
assert(strcmp(SKF2.label, "SKF2"))
assert(SKF2.n == n)
assert(SKF2.nu == nu)
assert(SKF2.ny == ny)
assert(SKF2.nf == size(SKF2.seq, 2))
assert(SKF2.nj == 2)
assert(isequal(SKF2.xkp1_est, x0))
assert(isequal(SKF2.Pkp1, P0))
assert(isequaln(SKF2.xk_est, nan(1, 1)))
assert(isequaln(SKF2.Pk, nan(1)))
assert(isequal(SKF2.rk, 1))

% Transition probabilities
epsilon = 0.05;
T = [1-epsilon epsilon; epsilon 1-epsilon];
assert(all(sum(T, 2) == 1))

% System indicator sequences
seq1 = {
    ones(1, nT+1);
    [ones(1, 20) 2*ones(1, nT+1-20)];  % equal to Gamma'
    [ones(1, 40) 2*ones(1, nT+1-40)];
    2*ones(1, nT+1);
 };
assert(isequal(seq1{2}, Gamma' + 1))

% Define MKF observer 1
seq = seq1;
nh = numel(seq);
%P0j = repmat({P0}, nh, 1);

% First, define with no initial state specified (should be set to zero)
% TODO: Allow independent P0 to be specified for each filter.
MKF1 = MKFObserverS(models,P0,seq,T,'MKF1');

assert(strcmp(MKF1.type, "MKF_S"))
assert(isequal(MKF1.models, models))
assert(isequal(MKF1.P0, P0))
assert(isequal(MKF1.seq, seq))
assert(isequal(MKF1.T, T))
assert(strcmp(MKF1.label, "MKF1"))
assert(MKF1.nh == nh)
assert(isequaln(MKF1.i, 0))
assert(isequal(MKF1.i_next, 1))
assert(MKF1.n == n)
assert(MKF1.nu == nu)
assert(MKF1.ny == ny)
assert(MKF1.nf == size(MKF1.seq{1}, 2))
assert(MKF1.nj == 2)
assert(isequal(MKF1.T, T))
assert(isequal(MKF1.xkp1_est, zeros(n, 1)))
assert(isequal(MKF1.Pkp1, P0))
assert(isequal(MKF1.r0, [1 1 1 2]'))
assert(isequal(MKF1.p_seq_g_Yk_init, ones(nh, 1) ./ nh))
assert(isequal(MKF1.rk, MKF1.r0))
assert(isequaln(MKF1.p_yk_g_seq_Ykm1, nan(nh, 1)))
assert(isequaln(MKF1.p_rk_g_rkm1, nan(nh, 1)))
assert(isequaln(MKF1.p_seq_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF1.xk_est, nan(n, 1)))
assert(isequaln(MKF1.Pk, nan(n, n)))
assert(isequaln(MKF1.yk_est, nan(ny, 1)))

% First, define with no initial state specified (should be set to zero)
% TODO: Allow independent P0 to be specified for each filter.
MKF2 = MKFObserverS(models,P0,seq,T,'MKF2');

assert(strcmp(MKF2.type, "MKF_S"))
assert(isequal(MKF2.models, models))
assert(isequal(MKF2.P0, P0))
assert(isequal(MKF2.seq, seq))
assert(isequal(MKF2.T, T))
assert(strcmp(MKF2.label, "MKF2"))
assert(MKF2.nh == nh)
assert(isequaln(MKF2.i, 0))
assert(isequal(MKF2.i_next, 1))
assert(MKF2.n == n)
assert(MKF2.nu == nu)
assert(MKF2.ny == ny)
assert(MKF2.nf == size(MKF2.seq{1}, 2))
assert(MKF2.nj == 2)
assert(isequal(MKF2.T, T))
assert(all(isnan(MKF2.xk_est)))
assert(isequaln(MKF2.Pk, nan))
assert(all(isnan(MKF2.yk_est)))
assert(isequal(MKF2.xkp1_est, zeros(n, 1)))
assert(isequal(MKF2.Pkp1, P0))
assert(isequal(MKF2.r0, [1 1 1 2]'))
assert(isequal(MKF2.p_seq_g_Yk_init, ones(nh, 1) ./ nh))
assert(isequal(MKF2.rk, MKF1.r0))
assert(isequaln(MKF2.p_yk_g_seq_Ykm1, nan(nh, 1)))
assert(isequaln(MKF2.p_rk_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF2.p_seq_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF2.xk_est, nan(n, 1)))
assert(isequaln(MKF2.Pk, nan(n, n)))
assert(isequaln(MKF2.yk_est, nan(ny, 1)))

% Redefine this time with initial conditions
MKF2 = MKFObserverS(models,P0,seq,T,'MKF2',x0);
assert(all(isnan(MKF2.xk_est)))
assert(all(isnan(MKF2.yk_est)))
assert(isequal(MKF2.xkp1_est, x0))
assert(isequal(MKF2.p_seq_g_Yk_init, ones(nh, 1) ./ nh))

% With initial prior shock values and probabilities
r0 = [2 2 2 2]';
MKF2 = MKFObserverS(models,P0,seq,T,'MKF2',x0,r0);
assert(isequal(MKF2.xkp1_est, x0))
assert(isequal(MKF2.rk, r0))
r0 = [1 1 1 2]';
p_seq_g_Yk_init = [0.6; 0.4];
MKF2 = MKFObserverS(A,B,C,Ts,P0,Q,R,seq,T,'MKF2',x0, ...
    r0,p_seq_g_Yk_init);
assert(isequal(MKF2.xkp1_est, x0))
assert(isequal(MKF2.rk, r0))
assert(isequal(MKF2.p_seq_g_Yk_init, p_seq_g_Yk_init))

% With default initial conditions
MKF2 = MKFObserverF(models,P0,seq,T,'MKF2');

% Choose observers to include in simulation
observers = {KF1, KF2, MKF1, MKF2, SKF1, SKF2};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Identify which observer to log MKF data for
f_mkf = 3;

% Simulate observers - without measurement noise (Y)
[Xk_est,Yk_est,DiagPk,Xkp1_est,Ykp1_est,DiagPkp1,MKF_K_obs,MKF_trP_obs, ...
    MKF_i,MKF_p_seq_g_Yk] = run_simulation_obs(Y,U,observers,f_mkf);

% Move prior estimates to correct time instants
Xkp1_est = [nan(1,n*n_obs); Xkp1_est(1:end-1,:)];
Ykp1_est = [nan(1,ny*n_obs); Ykp1_est(1:end-1,:)];

% Output estimation errors
E_obs = Y - Ykp1_est;

% Combine and display results
%sim_results1 = table(t,Gamma,U,X,Y,Ym,Xk_est,Yk_est,E_obs)

%figure(2); clf
%plot_obs_estimates(t,X,Xk_est,Y,Yk_est,obs_labels)

% Check KF1 was accurate before system switched
assert(max(abs(E_obs(t < 10, 1))) < 1e-5)

% Check MKF and SKF match KF1 before system switched
KF1_xkp1_est = Xkp1_est(t == 9.5, 1);
assert(isequal(abs(Xkp1_est(t == 9.5, :) - KF1_xkp1_est) < 0.0001, ...
    [true false true true true true]))
KF1_diagPk = sum(DiagPk(t == 9.5, 1));
KF1_diagPkp1 = sum(DiagPkp1(t == 9.5, 1));
assert(isequal(abs(DiagPk(t == 9.5, :) - KF1_diagPk) < 0.0001, ...
    [true false false true false true]))
assert(isequal(abs(DiagPkp1(t == 9.5, :) - KF1_diagPkp1) < 0.0001, ...
    [true false true false true true]))

% Check KF2 was accurate after system switched
assert(max(E_obs(t > 15, 2).^2) < 1e-3)

% Check MKF and SKF match KF2 after system switched
KF2_xk_est = Xk_est(t == 30, 2);
assert(isequal(abs(Xk_est(t == 30, :) - KF2_xk_est) < 0.0001, ...
    [false true false true false true]))
KF2_xkp1_est = Xk_est(t == 30, 2);
assert(isequal(abs(Xk_est(t == 30, :) - KF2_xkp1_est) < 0.0001, ...
    [false true false true false true]))
KF2_diagPk = sum(DiagPk(t == 30, 2));
assert(isequal(abs(DiagPk(t == 30, :) - KF2_diagPk) < 0.0001, ...
    [false true false true false true]))
KF2_diagPkp1 = sum(DiagPkp1(t == 30, 2));
% TODO: This is wrong
assert(isequal(abs(DiagPkp1(t == 30, :) - KF2_diagPkp1) < 0.0001, ...
    [false true true false true true]))

% Compute mean-squared error
mses = nanmean(E_obs.^2);

% Check MKF and SKF observer estimation errors
assert(isequal(round(mses, 4), [5.1728 0.4313 0.1296 0.1296 0.0660 0.0660]))

% Reset observer states to original initial conditions
KF1.reset()
KF2.reset()
MKF1.reset()
MKF2.reset()
SKF1.reset();
SKF2.reset();

assert(isequal(MKF1.P0, P0))
assert(isequal(MKF1.Pk, P0))
assert(isequal(MKF1.seq, seq))
assert(isequaln(MKF1.i, 0))
assert(isequal(MKF1.i_next, 1))
assert(isequal(MKF1.xkp1_est, zeros(n, 1)))
assert(MKF1.ykp1_est == 0)
assert(isequal(MKF1.rk, zeros(nh, 1)))
assert(isequaln(MKF1.p_yk_g_seq_Ykm1, nan(nh, 1)))
assert(isequaln(MKF1.p_rk_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF1.p_rk, nan(nh, 1)))
assert(isequaln(MKF1.p_seq_g_Ykm1, nan(nh, 1)))

assert(isequal(MKF2.P0, P0))
assert(isequal(MKF2.seq, seq))
assert(isequaln(MKF2.i, 0))
assert(isequal(MKF2.i_next, 1))
assert(isequaln(MKF2.xk_est, nan(n, 1)))
assert(isequaln(MKF2.Pk, nan(1)))
assert(isequaln(MKF2.yk_est, nan))
assert(isequal(MKF2.xkp1_est, zeros(n, 1)))
assert(isequal(MKF2.Pkp1, P0))
assert(MKF2.ykp1_est == 0)
assert(isequal(MKF2.rk, zeros(nh, 1)))
assert(isequaln(MKF2.p_yk_g_seq_Ykm1, nan(nh, 1)))
assert(isequaln(MKF2.p_rk_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF2.p_rk, nan(nh, 1)))
assert(isequaln(MKF2.p_seq_g_Ykm1, nan(nh, 1)))

% Redefine a new observer (identical to above)
MKF2_new = MKFObserverF(A,B,C,Ts,P0,Q,R,seq,T,'MKF2');
assert(isequaln(MKF2_new, MKF2))
MKF2_new.label = "MKF2_new";

% Make a copy
MKF2_copy = MKF2_new.copy();
assert(isequaln(MKF2_copy, MKF2_new))
MKF2_copy.label = "MKF2_copy";

% Choose observers to include in simulation
observers = {KF1, KF2, SKF2, MKF2, MKF2_new, MKF2_copy};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Identify which observer to log MKF data for
f_mkf = 4;

% Simulate observers - with measurement noise (Ym)
[Xk_est,Yk_est,DiagPk,Xkp1_est,Ykp1_est,DiagPkp1,MKF_K_obs,MKF_trP_obs, ...
    MKF_i,MKF_p_seq_g_Yk] = run_simulation_obs(Ym,U,observers,f_mkf);

% Move prediction estimates to correct time instants
Xk_km1_est = [nan(1,n*n_obs); Xkp1_est(1:end-1,:)];
Yk_km1_est = [nan(1,ny*n_obs); Ykp1_est(1:end-1,:)];

% Output estimation errors
E_obs = Y - Yk_est;

% Combine and display results
sim_results = table(t,Gamma,U,X,Y,Ym,Xk_est,Yk_est,E_obs);
writetable(sim_results, "results/test_MKFO_sim_results.csv");

% Display results from MKF observer
sim_results_MKF = [ ...
    table(t) ... 
    table(MKF_K_obs) ...
    table(MKF_trP_obs) ...
    table(MKF_i) ...
    table(MKF_p_seq_g_Yk) ...
];
writetable(sim_results_MKF, "results/test_MKFO_sim_results_MKF.csv");

% figure(3); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,obs_labels)

% % Plot difference between Xk_est and Xkp1_est for KF1
% figure(4); clf
% plot(t,X,"k--",t,Xk_km1_est(:,1),t,Xk_est(:,1),"Linewidth",2)
% grid on
% title("KF1")
% legend("x(k)","x(k|k-1)","x(k|k)",'Location','best')
% 
% % Plot difference between Xk_est and Xk_km1_est for KF2
% figure(5); clf
% plot(t,X,"k--",t,Xk_km1_est(:,2),t,Xk_est(:,2),"Linewidth",2)
% grid on
% title("KF2")
% legend("x(k)","x(k|k-1)","x(k|k)",'Location','best')

% Check final state estimates
test_Xk_est = [-1.191082  9.901224  9.901227  9.901227  9.915906  9.919637];
assert(isequal(round(Xk_est(t == t(end), :), 6), test_Xk_est))
% TODO: Why do the copies not produce identical simulation results?
% (see plot figure).

% Check final error covariance estimates
test_DiagPk = [0.015011  0.022516  0.022516  0.022516  0.022516  0.022516];
assert(isequal(round(DiagPk(t == t(end), :), 6), test_DiagPk))

% Compute mean-squared error
mses = nanmean(E_obs.^2);
%array2table(mses,'VariableNames',obs_labels)

% Check MKF observer estimation error
assert(round(mses(f_mkf), 4) == 0.0461)

% Check all observer estimation errors
assert(isequal(round(mses, 4), ...
    [3.8062 0.2694 0.0000 0.0461 0.0521 0.0126]))

% % Plot selected observers
% figure(4); clf
% plot_obs_estimates(t,X,X_est(:,[3 4]),Y,Y_est(:,[3 4]),obs_labels([3 4]))


%% Simulation test on 2x2 system

% Sample time
Ts = 1;

% NOTE: this is a previous version of the system with lower
% coupling (-0.2) and epsilon = [0.01; 0.01].

% Discrete time state space model
A = [ 0.8890       0     1 -0.2;
           0  0.8890  -0.2    1;
           0       0     1    0;
           0       0     0    1];
B = [    1 -0.2  0  0;
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

% Designate measured input and output signals
u_meas = [true; true; false; false];
y_meas = [true; true];

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);
nw = sum(~u_meas);

% Disturbance input (used by SKF observer)
Bw = B(:, ~u_meas);
nw = sum(~u_meas);

% RODD random variable parameters
epsilon = [0.01; 0.01];
sigma_M = [0.1; 0.1];
sigma_wp = [0.01 1; 0.01 1];

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

% Custom MKF test observers
% Devise a custom multi-model filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% this test simulation (t = t_shock)
% Multiple model filter 1
A2 = repmat({A}, 1, 3);
Bu2 = repmat({Bu}, 1, 3);
C2 = repmat({C}, 1, 3);
Du2 = repmat({Du}, 1, 3);
P0 = 1000*eye(n);
%P0_init = repmat({P0}, 1, 3);
Q2 = {diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,1)^2]), ...
      diag([0.01 0.01 sigma_wp(1,2)^2 sigma_wp(2,1)^2]), ...
      diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,2)^2])};
R2 = {diag(sigma_M.^2), diag(sigma_M.^2), diag(sigma_M.^2)};
seq = {zeros(1, nT+1); zeros(1, nT+1); zeros(1, nT+1); zeros(1, nT+1)};
seq{2}(t == t_shock(1)) = 1;  % shock 1
seq{3}(t == t_shock(2)) = 2;  % shock 2
seq{4}(t == t_shock(1)) = 1;  % both
seq{4}(t == t_shock(2)) = 2;
p_gamma = [1-epsilon epsilon]';
Z = [0 0; 1 0; 0 1; 1 1];  % combinations
p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
p_gamma = p_gamma ./ sum(p_gamma);  % normalized
T = repmat(p_gamma', 4, 1);
MKF3 = MKFObserver(A2,Bu2,C2,Ts,P0,Q2,R2,seq,T,'MKF3');
assert(MKF3.nh == 4)
MKF3F = MKFObserverF(A2,Bu2,C2,Ts,P0,Q2,R2,seq,T,'MKF3F');

seq = {zeros(1, nT+1)};
seq{1}(t == t_shock(1)) = 1;
seq{1}(t == t_shock(2)) = 2;
MKF4 = MKFObserver(A2,Bu2,C2,Ts,P0,Q2,R2,seq,T,'MKF4');
assert(MKF4.nh == 1)
MKF4F = MKFObserverF(A2,Bu2,C2,Ts,P0,Q2,R2,seq,T,'MKF4F');

% Define scheduled Kalman filter
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
% seq = sum(alpha .* 2.^(1:-1:0), 2)';
SKF = MKFObserverSched(A2,Bu2,C2,Ts,P0,Q2,R2,seq{1},"SKF");
SKFF = MKFObserverSchedF(A2,Bu2,C2,Ts,P0,Q2,R2,seq{1},"SKFF");

% Choose observers to test
observers = {MKF3, MKF4, SKF, MKF3F, MKF4F, SKFF};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Simulate system
X = zeros(nT+1,n);
Y = zeros(nT+1,ny);
xk = zeros(n,1);

for i = 1:nT+1

    % Inputs
    uk = U_sim(i,:)';

    % Compute y(k)
    yk = C*xk + D*uk;

    % Store results
    X(i, :) = xk';
    Y(i, :) = yk';

    % Compute x(k+1)
    xk = A*xk + B*uk;

end

% Check simulation output is correct
[Y2, t, X2] = lsim(Gpss, U_sim, t);
assert(isequal(X, X2))
assert(isequal(Y, Y2))

% Choose measurement noise for plant
sigma_MP = [0; 0];  % Set to zero for testing
Y_m = Y + sigma_MP'.*randn(nT+1, ny);

% Identify which observer to log MKF data for
f_mkf = 1;

% Simulate observers
[Xk_est,Yk_est,DiagPk,Xkp1_est,Ykp1_est,DiagPkp1,MKF_K_obs,MKF_trP_obs, ...
    MKF_i,MKF_p_seq_g_Yk] = run_simulation_obs(Y,U,observers,f_mkf);

% Move prediction estimates to correct time instants
Xk_km1_est = [nan(1,n*n_obs); Xkp1_est(1:end-1,:)];
Yk_km1_est = [nan(1,ny*n_obs); Ykp1_est(1:end-1,:)];

% Output estimation errors
E_obs_ykm1 = repmat(Y,1,n_obs) - Yk_km1_est;
E_obs_yk = repmat(Y,1,n_obs) - Yk_est;

% Plot observer estimates
% figure(7); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,obs_labels)

% Check final state predition estimates
test_X_est = [-1.801802  9.009008  1.000000  1.000000 -1.801802 ...
    9.009008  1.000000  1.000000 -1.801802  9.009008  1.000000  1.000000];
assert(isequal(round(Xkp1_est(t == t(end), :), 6), [test_X_est test_X_est]))

% Check final error covariance estimates - P(k+1|k)
% TODO: Haven't checked if these are correct.
test_DiagPkp1 = [ 0.092947  0.092947  0.002086  0.002086  0.092947 ...
    0.092947  0.002086  0.002086  0.092947  0.092947  0.002086  0.002086];
assert(isequal(round(DiagPkp1(t == t(end), 1:12), 6), test_DiagPkp1))

% Check final error covariance estimates - P(k+1|k)
% TODO: Haven't checked if these are correct.
test_DiagPkp1 = [ 0.083317  0.083317  0.001986  0.001986  0.083317 ...
    0.083317  0.001986  0.001986  0.092947  0.092947  0.002086  0.002086];
assert(isequal(round(DiagPkp1(t == t(end), 13:24), 6), test_DiagPkp1))

% Check final error covariance estimates - P(k|k)
% TODO: Haven't checked if these are correct.
test_DiagPk = [ 0.083317  0.083317  0.001986  0.001986  0.083317 ...
    0.083317  0.001986  0.001986  0.083317  0.083317  0.001986  0.001986];
assert(isequal(round(DiagPk(t == t(end), 13:24), 6), test_DiagPk))

% Compute mean-squared error (use prediction estimates y(k|k-1)
% because all observers compute these)
mses_y12 = nanmean(E_obs_ykm1.^2);
mses_yi = array2table(nanmean(reshape(mses_y12,ny,n_obs)), ...
    'VariableNames',obs_labels);

% Check observer estimation errors
assert(isequal(round(mses_yi.Variables, 6), ...
    [0.000817 0.000128 0.000128 0.000817 0.000128 0.000128]))
% Should be the same for both sets of observers

% Compute mean-squared error with updated estimates y(k|k)
mses_y12 = nanmean(E_obs_yk(:, 7:12).^2);
mses_yi = array2table(nanmean(reshape(mses_y12,ny,3)), ...
    'VariableNames',obs_labels(4:6));

% Check observer estimation errors
assert(isequal(round(mses_yi.Variables, 7), ...
    [0.0003349  0.0000196  0.0000196]))
% Should be the same for MKF4F and SKFF

% % Display results of last simulation
% 
% X_est = sim_results.X_est;
% E_obs = sim_results.E_obs;
% K_obs = sim_results.K_obs;
% trP_obs = sim_results.trP_obs;
% 
% table(t,alpha,U,Wp,X,Y,Y_m,X_est,Y_est,E_obs)
% 
% % Display gains and trace of covariance matrix
% table(t, cell2mat(K_obs), cell2mat(trP_obs), ...
%     'VariableNames', {'t', 'K{1}, K{2}', 'trace(P{1}), trace(P{2})'})
% 
% % Show table of mean-squared errors
% table(MSE.keys', cell2mat(MSE.values'), ...
%     'VariableNames', {'Observer', 'MSE'})



%% Test copy methods

% Define system

% Sample period
Ts = 0.5;

% Discrete time state space models
% Model #1
A1 = 0.7;
B1 = 1;
C1 = 0.3;
D1 = 0;
Gpss1 = ss(A1,B1,C1,D1,Ts);

% Model #2
A2 = 0.9;
B2 = 1;
C2 = -0.3;  % -ve gain!
D2 = 0;
Gpss2 = ss(A2,B2,C2,D2,Ts);

% Dimensions
n = size(A1, 1);
nu = size(B1, 2);
ny = size(C1, 1);

% Check dimensions
assert(isequal(size(A1), size(A2)))
assert(isequal(size(B1), size(B2)))
assert(isequal(size(C1), size(C2)))
assert(isequal(size(D1), size(D2)))

% Define system models
A = {A1, A2};
B = {B1, B2};
C = {C1, C2};
D = {D1, D2};

% Observer parameters (same for all observers)
P0 = 10000;
x0 = 0.5;
Q1 = 0.01;
R1 = 0.1^2;
Q2 = 0.01;
R2 = 0.1^2;

% Switching parameters
Q = {Q1,Q2};
R = {R1,R2};

% Transition probabilities
epsilon = 0.05;
T = [1-epsilon epsilon; epsilon 1-epsilon];
assert(all(sum(T, 2) == 1))

% System indicator sequences
nT = 60;
seq1 = {
    zeros(1, nT+1);
    [zeros(1, 20) ones(1, nT+1-20)];  % equal to Gamma'
    [zeros(1, 40) ones(1, nT+1-40)];
    ones(1, nT+1);
 };

% Define MKF observer
seq = seq1;
nh = numel(seq);
%P0j = repmat({P0}, nh, 1);

% Define multi-model observer with initial conditions
MKF = MKFObserverF(A,B,C,Ts,P0,Q,R,seq,T,'MKF',x0);

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
assert(isequaln(MKF_copy.filters{1}, MKF.filters{1}))
assert(isequaln(MKF_copy.filters{2}, MKF.filters{2}))

% Check deep copy was made
% TODO: This is not working
%assert(MKF_copy.filters{1} ~= MKF.filters{1})  % must not be same object
%assert(MKF_copy.filters{2} ~= MKF.filters{2})

MKF.label = "New name";
assert(~isequal(MKF_copy.label, "New name"))

MKF.filters{1}.x0 = 99;
%assert(~isequal(MKF_copy.filters{1}.x0, 99))

%END


% function [X, Y, Ym] = run_simulation_sys(A,B,C,D,U,V,Gamma,nT)
% % Simulate switching system
% 
%     [n, ~, ny] = check_dimensions(A{1}, B{1}, C{1}, D{1});
%     X = zeros(nT+1,n);
%     Y = zeros(nT+1,ny);
%     Ym = zeros(nT+1,ny);
%     xk = zeros(n,1);
% 
%     for i=1:nT+1
% 
%         % Switch system
%         j = Gamma(i) + 1;
% 
%         % Inputs
%         uk = U(i,:)';
% 
%         % Compute y(k)
%         yk = C{j}*xk + D{j}*uk;
%         yk_m = yk + V(i);
% 
%         % Store results
%         X(i, :) = xk';
%         Y(i, :) = yk';
%         Ym(i, :) = yk_m';
% 
%         % Compute x(k+1)
%         xk = A{j}*xk + B{j}*uk;
% 
%     end
% end


function [Xk_est,Yk_est,DiagPk,Xkp1_est,Ykp1_est,DiagPkp1,MKF_K_obs,MKF_trP_obs, ...
    MKF_i,MKF_p_seq_g_Yk] = run_simulation_obs(Ym,U,observers,f_mkf)
% Simulate observers

    nT = size(Ym, 1) - 1;
    ny = size(Ym, 2);
    n_obs = numel(observers);
    n = size(observers{1}.xkp1_est, 1);

    obs_mkf = observers{f_mkf};
    nhers = size(obs_mkf.seq, 1);

    Xk_est = nan(nT+1, n*n_obs);
    Yk_est = nan(nT+1, ny*n_obs);
    Xkp1_est = nan(nT+1, n*n_obs);
    Ykp1_est = nan(nT+1, ny*n_obs);
    DiagPk = nan(nT+1, n*n_obs);
    DiagPkp1 = nan(nT+1, n*n_obs);
    MKF_K_obs = cell(nT+1, n*nhers);
    MKF_trP_obs = nan(nT+1, nhers);
    MKF_i = nan(nT+1, 2);
    MKF_p_seq_g_Yk = nan(nT+1, nhers);

    for i = 1:nT+1

        yk = Ym(i, :)';
        uk = U(i, :)';

        % Update observers
        for f = 1:n_obs
            obs = observers{f};
            obs.update(yk, uk);
            if f == f_mkf
                MKF_i(i, :) = obs.i;
                MKF_p_seq_g_Yk(i, :) = obs.p_seq_g_Yk';
                for j = 1:obs.nh
                    switch obs.type
                        case {"KF", "SKF", "MKF"}
                            MKF_K_obs{i, j} = obs.filters{j}.K';
                            MKF_trP_obs(i, j) = trace(obs.filters{j}.Pkp1);
                        case {"KFF", "SKFF", "MKFF"}
                            MKF_K_obs{i, j} = obs.filters{j}.Kf';
                            MKF_trP_obs(i, j) = trace(obs.filters{j}.Pkp1);
                    end
                end
            end
            if isprop(obs,'xk_est')
                xk_est(1, (f-1)*n+1:f*n) = obs.xk_est';
                yk_est(1, (f-1)*ny+1:f*ny) = obs.yk_est';
            end
            if isprop(obs,'xkp1_est')
                xkp1_est(1, (f-1)*n+1:f*n) = obs.xkp1_est';
                ykp1_est(1, (f-1)*ny+1:f*ny) = obs.ykp1_est';
            end
            if isprop(obs,'P')
                diagPkp1(1, (f-1)*n+1:f*n) = diag(obs.P)';
            else
                if isprop(obs,'Pk')
                    diagPk(1, (f-1)*n+1:f*n) = diag(obs.Pk)';
                end
                diagPkp1(1, (f-1)*n+1:f*n) = diag(obs.Pkp1)';
            end
        end

        % Record observer estimates
        Xk_est(i, :) = xk_est;
        Yk_est(i, :) = yk_est;
        Xkp1_est(i, :) = xkp1_est;
        Ykp1_est(i, :) = ykp1_est;
        DiagPk(i, :) = diagPk;
        DiagPkp1(i, :) = diagPkp1;

    end
end
