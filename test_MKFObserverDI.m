% Test classes MKFObserverDI and MKFObserverSched
% TODO: At the moment this just replicates tests in MKFObserver.m
% Need to setup additional tests to check when d > 1.

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
A = {A1, A2};
B = {B1, B2};
C = {C1, C2};
D = {D1, D2};

% Input disturbance variance
%sigma_w = 0.1;
sigma_w = 0;

% Process noise std. dev.
sigma_W = [0; 0];

% Measurement noise std. dev.
sigma_M = 0.1;

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
[X, Y, Ym] = run_simulation_sys(A,B,C,D,U,V,Gamma,nT);

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
Q1 = 0.01;
R1 = 0.1^2;
Q2 = 0.01;
R2 = 0.1^2;
assert(isequal(size(Q1), size(Q2)))
assert(isequal(size(R1), size(R2)))

% Standard Kalman filters
KF1 = KalmanFilter(A1,B1,C1,Ts,P0,Q1,R1,'KF1',x0);
KF2 = KalmanFilter(A2,B2,C2,Ts,P0,Q2,R2,'KF2',x0);

% Define observers with a switching system
Q = {Q1,Q2};
R = {R1,R2};

% Define scheduled MKF filter
seq = Gamma';
SKF = MKFObserverSched(A,B,C,Ts,P0,Q,R,seq,"SKF1",x0);

assert(strcmp(SKF.type, "SKF"))
assert(isequal(SKF.A, A))
assert(isequal(SKF.B, B))
assert(isequal(SKF.C, C))
assert(isequal(SKF.Ts, Ts))
assert(isequal(SKF.P0, P0))
assert(isequal(SKF.P, P0))
assert(isequal(SKF.Q, Q))
assert(isequal(SKF.R, R))
assert(isequal(SKF.seq, seq))
assert(strcmp(SKF.label, "SKF1"))
assert(SKF.n_filt == 1)
assert(isa(SKF.filter, 'KalmanFilter'))
assert(strcmp(SKF.filter.label, 'KF'))
assert(SKF.n == n)
assert(SKF.nu == nu)
assert(SKF.ny == ny)
assert(SKF.f == size(SKF.seq, 2))
assert(SKF.nj == 2)
assert(isequal(SKF.xkp1_est, x0))
assert(SKF.ykp1_est == C{1}*x0)
assert(isequal(SKF.gamma_k, 0))

% Transition probabilities
epsilon = 0.05;
T = [1-epsilon epsilon; epsilon 1-epsilon];
assert(all(sum(T, 2) == 1))

% System indicator sequences
seq1 = {
    zeros(1, nT+1);
    [zeros(1, 20) ones(1, nT+1-20)];  % equal to Gamma'
    [zeros(1, 40) ones(1, nT+1-40)];
    ones(1, nT+1);
 };
assert(isequal(seq1{2}, Gamma'))

% Define MKF observer 1
seq = seq1;
n_filt = numel(seq);
%P0j = repmat({P0}, n_filt, 1);
d = 1;

% First, define with no initial state specified (should be set to zero)
% TODO: Allow independent P0 to be specified for each filter.
MKF1 = MKFObserverDI(A,B,C,Ts,P0,Q,R,seq,T,d,'MKF1');

assert(strcmp(MKF1.type, "MKF_DI"))
assert(isequal(MKF1.A, A))
assert(isequal(MKF1.B, B))
assert(isequal(MKF1.C, C))
assert(isequal(MKF1.Ts, Ts))
assert(isequal(MKF1.P0, P0))
assert(isequal(MKF1.P, P0))
assert(isequal(MKF1.Q, Q))
assert(isequal(MKF1.R, R))
assert(isequal(MKF1.seq, seq))
assert(isequal(MKF1.T, T))
assert(isequal(MKF1.d, d))
assert(strcmp(MKF1.label, "MKF1"))
assert(MKF1.n_filt == n_filt)
assert(MKF1.n_filt == numel(MKF1.filters))
assert(strcmp(MKF1.filters{n_filt}.label, 'MKF14'))
assert(isequaln(MKF1.i, [0 0]))
assert(isequal(MKF1.i_next, int16([1 1])))
assert(MKF1.n == n)
assert(MKF1.nu == nu)
assert(MKF1.ny == ny)
assert(MKF1.f == size(MKF1.seq{1}, 2))
assert(MKF1.nj == 2)
assert(isequal(MKF1.T, T))
assert(isequal(MKF1.xkp1_est, zeros(n, 1)))
assert(MKF1.ykp1_est == 0)
assert(isequal(MKF1.gamma_k, zeros(n_filt, 1)))
assert(isequal(MKF1.p_yk_g_seq_Ykm1, zeros(n_filt, 1)))
assert(isequal(MKF1.p_gammak_g_Ykm1, zeros(n_filt, 1)))
assert(isequal(MKF1.p_gamma_k, [0.95 0.95 0.95 0.95]'))
assert(isequal(MKF1.p_seq_g_Ykm1, zeros(n_filt, 1)))

% Redefine this time with initial conditions
MKF1 = MKFObserverDI(A,B,C,Ts,P0,Q,R,seq,T,d,'MKF1',x0);
assert(isequal(MKF1.xkp1_est, x0))
assert(isequal(MKF1.ykp1_est, C{1} * x0))
gamma0 = 0;
MKF1 = MKFObserverDI(A,B,C,Ts,P0,Q,R,seq,T,d,'MKF1',x0,gamma0);
assert(isequal(MKF1.xkp1_est, x0))
assert(isequal(MKF1.ykp1_est, C{1} * x0))
assert(isequal(MKF1.gamma_k, zeros(n_filt, 1)))
gamma0 = zeros(n_filt, 1);
gamma0(end) = 1;
MKF1 = MKFObserverDI(A,B,C,Ts,P0,Q,R,seq,T,d,'MKF1',x0,gamma0);
assert(isequal(MKF1.xkp1_est, x0))
assert(isequal(MKF1.ykp1_est, C{1} * x0))
assert(isequal(MKF1.gamma_k, gamma0))

% With default initial conditions
MKF1 = MKFObserverDI(A,B,C,Ts,P0,Q,R,seq,T,d,'MKF1');

% Choose observers to include in simulation
observers = {KF1, KF2, MKF1, SKF};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Identify which observer to log MKF data for
f_mkf = 3;

% Simulate observers - without measurement noise (Y)
[Xkp1_est,Ykp1_est,DiagP,MKF_K_obs,MKF_trP_obs,MKF_i,MKF_p_seq_g_Yk] = ...
    run_simulation_obs(Y,U,observers,f_mkf);

% Move estimates to correct time instants
X_est = [nan(1,n*n_obs); Xkp1_est(1:end-1,:)];
Y_est = [nan(1,ny*n_obs); Ykp1_est(1:end-1,:)];
P = [nan(1,n*n_obs); DiagP(1:end-1,:)];

% Output estimation errors
E_obs = Y - Y_est;

% Combine and display results
%sim_results1 = table(t,Gamma,U,X,Y,Ym,X_est,Y_est,E_obs)

%figure(2); clf
%plot_obs_estimates(t,X,X_est,Y,Y_est,obs_labels)

% Check KF1 was accurate before system switched
assert(nanmean(E_obs(t < 10, 1).^2) < 0.0001)

% Check MKF and SKF match KF1 before system switched
KF1_x_est = X_est(t == 9.5, 1);
assert(isequal(abs(X_est(t == 9.5, :) - KF1_x_est) < 0.0001, ...
    [true false true true]))
KF1_diagP = sum(DiagP(t == 9.5, 1));
assert(isequal(abs(DiagP(t == 9.5, :) - KF1_diagP) < 0.0001, ...
    [true false true true]))

% Check KF2 was accurate after system switched
assert(mean(E_obs(t > 12, 2).^2) < 0.001)

% Check MKF and SKF match KF2 after system switched
KF2_x_est = X_est(t == 30, 2);
assert(isequal(abs(X_est(t == 30, :) - KF2_x_est) < 0.0001, ...
    [false true true true]))
KF2_diagP = sum(DiagP(t == 30, 2));
assert(isequal(abs(DiagP(t == 30, :) - KF2_diagP) < 0.0001, ...
    [false true true true]))

% Compute mean-squared error
mses = nanmean(E_obs.^2);

% Check MKF and SKF observer estimation errors
assert(isequal(round(mses, 4), [5.1728 0.4313 0.1296 0.0660]))

% Reset observer states to original initial conditions
KF1.reset()
KF2.reset()
MKF1.reset()
SKF.reset();

assert(isequal(MKF1.P0, P0))
assert(isequal(MKF1.P, P0))
assert(isequal(MKF1.seq, seq))
assert(isequaln(MKF1.i, [0 0]))
assert(isequal(MKF1.i_next, int16([1 1])))
assert(isequal(MKF1.xkp1_est, zeros(n, 1)))
assert(MKF1.ykp1_est == 0)
assert(isequal(MKF1.gamma_k, zeros(n_filt, 1)))
assert(isequal(MKF1.p_yk_g_seq_Ykm1, zeros(n_filt, 1)))
assert(isequal(MKF1.p_gammak_g_Ykm1, zeros(n_filt, 1)))
assert(isequal(MKF1.p_gamma_k, [0.95 0.95 0.95 0.95]'))
assert(isequal(MKF1.p_seq_g_Ykm1, zeros(n_filt, 1)))

% Redefine a new observer (identical to above)
MKF1_new = MKFObserverDI(A,B,C,Ts,P0,Q,R,seq,T,d,'MKF1');
assert(isequaln(MKF1_new, MKF1))
MKF1_new.label = "MKF1_new";

% Make a copy
MKF1_copy = MKF1_new.copy();
assert(isequaln(MKF1_copy, MKF1_new))
MKF1_copy.label = "MKF1_copy";

% Choose observers to include in simulation
observers = {KF1, KF2, SKF, MKF1, MKF1_new, MKF1_copy};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Identify which observer to log MKF data for
f_mkf = 4;

% Simulate observers - with measurement noise (Ym)
[Xkp1_est,Ykp1_est,DiagP,MKF_K_obs,MKF_trP_obs,MKF_i,MKF_p_seq_g_Yk] = ...
    run_simulation_obs(Ym,U,observers,f_mkf);

% Move estimates to correct time instants
X_est = [nan(1,n*n_obs); Xkp1_est(1:end-1,:)];
Y_est = [nan(1,ny*n_obs); Ykp1_est(1:end-1,:)];

% Output estimation errors
E_obs = Y - Y_est;

% Combine and display results
sim_results = table(t,Gamma,U,X,Y,Ym,X_est,Y_est,E_obs);
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
% plot_obs_estimates(t,X,X_est,Y,Y_est,obs_labels)

% Check final state estimates
test_X_est = [0.202984  9.838603  9.838607  9.838607  9.854841  9.807889];
assert(isequal(round(X_est(t == t(end), :), 6), test_X_est))
% TODO: Why do the copies not produce identical simulation results?
% (see plot figure).

% Check final error covariance estimates
test_DiagP = [0.017355  0.028238  0.028238  0.028238  0.028238  0.028238];
assert(isequal(round(DiagP(t == t(end), :), 6), test_DiagP))

% Compute mean-squared error
mses = nanmean(E_obs.^2);
%array2table(mses,'VariableNames',obs_labels)

% Check MKF observer estimation error
assert(round(mses(f_mkf), 4) == 0.1335)

% Check all observer estimation errors
assert(isequal(round(mses, 4), ...
    [5.1749  0.4074  0.0685  0.1335  0.1387  0.0877]))
% TODO: Why do the copies not produce identical simulation results?
% (see plot figure).

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
Z = [0 0; 1 0; 0 1];  % combinations
p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
p_gamma = p_gamma ./ sum(p_gamma);  % normalized
T = repmat(p_gamma', 3, 1);
d = 1;
MKF3 = MKFObserverDI(A2,Bu2,C2,Ts,P0,Q2,R2,seq,T,d,'MKF3');
assert(MKF3.n_filt == 4)

seq = {zeros(1, nT+1)};
seq{1}(t == t_shock(1)) = 1;
seq{1}(t == t_shock(2)) = 2;
MKF4 = MKFObserverDI(A2,Bu2,C2,Ts,P0,Q2,R2,seq,T,d,'MKF4');
assert(MKF4.n_filt == 1)

% Define scheduled Kalman filter
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
% seq = sum(alpha .* 2.^(1:-1:0), 2)';
SKF = MKFObserverSched(A2,Bu2,C2,Ts,P0,Q2,R2,seq{1},"SKF");

% Choose observers to test
observers = {MKF3, MKF4, SKF};
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
[Xkp1_est,Ykp1_est,DiagP,MKF_K_obs,MKF_trP_obs,MKF_i,MKF_p_seq_g_Yk] = ...
    run_simulation_obs(Y_m,U,observers,f_mkf);

% Plot observer estimates
% figure(5); clf
% plot_obs_estimates(t,X,Xkp1_est,Y,Ykp1_est,obs_labels)

% Check final state estimates
test_X_est = [-1.801802  9.009008  1.000000  1.000000 -1.801802 ...
    9.009008  1.000000  1.000000 -1.801802  9.009008  1.000000  1.000000];
assert(isequal(round(Xkp1_est(t == t(end), :), 6), test_X_est))

% Check final error covariance estimates
% TODO: Haven't checked if these are correct.
test_DiagP = [ 0.092947  0.092947  0.002086  0.002086  0.092947 ...
    0.092947  0.002086  0.002086  0.092947  0.092947  0.002086  0.002086];
assert(isequal(round(DiagP(t == t(end), :), 6), test_DiagP))

% Compute mean-squared errors
MSE = containers.Map();
for i = 1:n_obs
    MSE(observers{i}.label) = mean((Ykp1_est(:, i*ny-1:i*ny) - Y).^2);
end
%fprintf("%d, %s: %f\n", i, obs.label, mean((Y_est - Y).^2))

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

% Results on Nov 8 after reverting back the Bayesian updating
MSE_test_values = containers.Map(...
 {'MKF3',              'MKF4',              'SKF'             }, ...
 {[0.000707 0.000348], [0.000234 0.000083], [0.000234 0.000083]} ...
);

for label = MSE.keys
    %fprintf("%s: %f, %f (%f, %f)\n", label{1}, MSE(label{1}), MSE_test_values(label{1}))
    assert(isequal(round(MSE(label{1}), 6), MSE_test_values(label{1})))
end


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
n_filt = numel(seq);
%P0j = repmat({P0}, n_filt, 1);
d = 1;

% Define multi-model observer with initial conditions
MKF = MKFObserverDI(A,B,C,Ts,P0,Q,R,seq,T,d,'MKF',x0);

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


function [X, Y, Ym] = run_simulation_sys(A,B,C,D,U,V,Gamma,nT)
% Simulate switching system

    [n, ~, ny] = check_dimensions(A{1}, B{1}, C{1}, D{1});
    X = zeros(nT+1,n);
    Y = zeros(nT+1,ny);
    Ym = zeros(nT+1,ny);
    xk = zeros(n,1);

    for i=1:nT+1

        % Switch system
        j = Gamma(i) + 1;

        % Inputs
        uk = U(i,:)';

        % Compute y(k)
        yk = C{j}*xk + D{j}*uk;
        yk_m = yk + V(i);

        % Store results
        X(i, :) = xk';
        Y(i, :) = yk';
        Ym(i, :) = yk_m';

        % Compute x(k+1)
        xk = A{j}*xk + B{j}*uk;

    end
end

% TODO: Can't the function run_test_simulation.m be used here?
function [Xkp1_est,Ykp1_est,DiagP,MKF_K_obs,MKF_trP_obs,MKF_i,MKF_p_seq_g_Yk] = ...
    run_simulation_obs(Ym,U,observers,f_mkf)
% Simulate observers

    nT = size(Ym, 1) - 1;
    ny = size(Ym, 2);
    n_obs = numel(observers);
    n = size(observers{1}.xkp1_est, 1);

    obs_mkf = observers{f_mkf};
    n_filters = size(obs_mkf.seq, 1);

    Xkp1_est = zeros(nT+1, n*n_obs);
    Ykp1_est = zeros(nT+1, ny*n_obs);
    DiagP = zeros(nT+1, n*n_obs);
    MKF_K_obs = cell(nT+1, n*n_filters);
    MKF_trP_obs = nan(nT+1, n_filters);
    MKF_i = nan(nT+1, 2);
    MKF_p_seq_g_Yk = nan(nT+1, n_filters);

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
                for j = 1:obs.n_filt
                    MKF_K_obs{i, j} = obs.filters{j}.K';
                    MKF_trP_obs(i, j) = trace(obs.filters{j}.P);
                end
            end
            xkp1_est(1, (f-1)*n+1:f*n) = obs.xkp1_est';
            ykp1_est(1, (f-1)*ny+1:f*ny) = obs.ykp1_est';
            diagP(1, (f-1)*n+1:f*n) = diag(obs.P)';
        end

        % Record observer estimates
        Xkp1_est(i, :) = xkp1_est;
        Ykp1_est(i, :) = ykp1_est;
        DiagP(i, :) = diagP;

    end
end
