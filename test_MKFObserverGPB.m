% Test the following multi-model observer classes:
%  - MKFObserverAMM, MKFObserverGPB1, and MKFObserverGPB2

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
%Gpss1 = ss(A1,B1,C1,D1,Ts);

% Model #2
A2 = 0.9;
B2 = 1;
C2 = -0.3;  % -ve gain!
D2 = 0;
%Gpss2 = ss(A2,B2,C2,D2,Ts);

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

% Actual system switching sequence
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
x0 = 0.5;  % optional
y0 = 0.1;  % optional
Q1 = 0.01;
R1 = 0.1^2;
Q2 = 0.01;
R2 = 0.1^2;
assert(isequal(size(Q1), size(Q2)))
assert(isequal(size(R1), size(R2)))

% Standard Kalman filters
KF1 = KalmanFilterF(A1,B1,C1,Ts,P0,Q1,R1,'KF1',x0);
KF2 = KalmanFilterF(A2,B2,C2,Ts,P0,Q2,R2,'KF2',x0);

% Define observers with a switching system
Q = {Q1,Q2};
R = {R1,R2};

% Define scheduled MKF filter
seq = Gamma';
SKF = MKFObserverSchedF(A,B,C,Ts,P0,Q,R,seq,"SKF1",x0);

assert(strcmp(SKF.type, "SKFF"))
assert(isequal(SKF.A, A))
assert(isequal(SKF.B, B))
assert(isequal(SKF.C, C))
assert(isequal(SKF.Ts, Ts))
assert(isequal(SKF.P0, P0))
assert(isequaln(SKF.Pk, nan))
assert(isequal(SKF.Pkp1, P0))
assert(isequal(SKF.Q, Q))
assert(isequal(SKF.R, R))
assert(isequal(SKF.seq, seq))
assert(strcmp(SKF.label, "SKF1"))
assert(SKF.n_filt == 1)
assert(isa(SKF.filter, 'KalmanFilterF'))
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
T = [0.95 0.05; 0.01 0.99];
assert(all(sum(T, 2) == 1))

% System indicator sequences
seq1 = {
    zeros(1, nT+1);
    [zeros(1, 20) ones(1, nT+1-20)];  % equal to Gamma'
    [zeros(1, 40) ones(1, nT+1-40)];
    ones(1, nT+1);
 };
assert(isequal(seq1{2}, Gamma'))

% Define MKF observers

% Two system models
m1.A = A1;
m1.B = B1;
m1.C = C1;
m1.Q = Q1;
m1.R = R1;
m2.A = A2;
m2.B = B2;
m2.C = C2;
m2.Q = Q2;
m2.R = R2;
models = {m1, m2};

% 1. AMM
MKF_AMM1 = MKFObserverAMM(models,Ts,P0,"MKF_AMM1");
% MKF_AMM1 = MKFObserverAMM(models,Ts,P0,"MKF_AMM1",x0,y0,p_seq_g_Yk_init);

assert(strcmp(MKF_AMM1.type, "MKF_AMM"))
assert(isequal(MKF_AMM1.models, models))
assert(isequal(MKF_AMM1.Ts, Ts))
assert(isequal(MKF_AMM1.P0, P0))
assert(strcmp(MKF_AMM1.label, "MKF_AMM1"))
assert(MKF_AMM1.n_filt == 2)
assert(MKF_AMM1.n == n)
assert(MKF_AMM1.nu == nu)
assert(MKF_AMM1.ny == ny)
assert(MKF_AMM1.nj == 2)
assert(isequal(MKF_AMM1.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_AMM1.Pkp1, P0))
assert(isequal(MKF_AMM1.ykp1_est, 0))
assert(isequaln(MKF_AMM1.xk_est, nan(n, 1)))
assert(isequaln(MKF_AMM1.Pk, nan(n)))
assert(isequaln(MKF_AMM1.yk_est, nan(ny, 1)))
assert(isequaln(MKF_AMM1.x0, zeros(n, 1)))
assert(isequaln(MKF_AMM1.y0, zeros(ny, 1)))
assert(isequal(MKF_AMM1.p_seq_g_Yk_init, [0.5 0.5]'))

% Redefine this time with initial conditions
MKF_AMM1x0 = MKFObserverAMM(models,Ts,P0,'MKF_AMM1x0',x0);
assert(isequaln(MKF_AMM1x0.x0, x0))
assert(isequal(MKF_AMM1x0.xkp1_est, x0))
assert(isequal(MKF_AMM1x0.Pkp1, P0))
assert(isequal(MKF_AMM1x0.ykp1_est, C1 * x0))
assert(isequal(MKF_AMM1x0.p_seq_g_Yk_init, [0.5; 0.5]))
MKF_AMM1x0y0 = MKFObserverAMM(models,Ts,P0,'MKF_AMM1x0',x0,y0);
assert(isequaln(MKF_AMM1x0y0.x0, x0))
assert(isequaln(MKF_AMM1x0y0.y0, y0))
assert(isequal(MKF_AMM1x0y0.xkp1_est, x0))
assert(isequal(MKF_AMM1x0y0.ykp1_est, y0))

% Also with initial prior probability values
p_seq_g_Yk_init = [0.6; 0.4];
MKF_AMM1x0y0p0 = MKFObserverAMM(models,Ts,P0,"MKF_GPB1x0",x0,y0, ...
    p_seq_g_Yk_init);
assert(isequaln(MKF_AMM1x0y0p0.x0, x0))
assert(isequaln(MKF_AMM1x0y0p0.y0, y0))
assert(isequaln(MKF_AMM1x0y0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))
assert(isequal(MKF_AMM1x0y0p0.xkp1_est, x0))
assert(isequal(MKF_AMM1x0y0p0.ykp1_est, y0))
assert(isequal(MKF_AMM1x0y0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))

% % 2. GPB1
MKF_GPB1 = MKFObserverGPB1(models,Ts,P0,T,"MKF_GPB1");

assert(strcmp(MKF_GPB1.type, "MKF_GPB1"))
assert(isequal(MKF_GPB1.models, models))
assert(isequal(MKF_GPB1.Ts, Ts))
assert(isequal(MKF_GPB1.P0, P0))
assert(isequal(MKF_GPB1.T, T))
assert(strcmp(MKF_GPB1.label, "MKF_GPB1"))
assert(MKF_GPB1.n_filt == 2)
assert(MKF_GPB1.n == n)
assert(MKF_GPB1.nu == nu)
assert(MKF_GPB1.ny == ny)
assert(MKF_GPB1.nj == 2)
assert(isequal(MKF_GPB1.T, T))
assert(isequal(MKF_GPB1.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_GPB1.Pkp1, P0))
assert(isequal(MKF_GPB1.ykp1_est, 0))
assert(isequaln(MKF_GPB1.xk_est, nan(n, 1)))
assert(isequaln(MKF_GPB1.Pk, nan(n)))
assert(isequaln(MKF_GPB1.yk_est, nan(ny, 1)))
assert(isequaln(MKF_GPB1.gamma_k, [0 1]'))
assert(isequaln(MKF_GPB1.p_yk_g_seq_Ykm1, nan(MKF_GPB1.n_filt, 1)))
assert(isequaln(MKF_GPB1.p_gammak_g_Ykm1, nan(MKF_GPB1.n_filt, 1)))
assert(isequaln(MKF_GPB1.p_gamma_k_g_gamma_km1, [0.95 0.99]'))  % This doesn't make sense
assert(isequaln(MKF_GPB1.p_seq_g_Ykm1, nan(MKF_GPB1.n_filt, 1)))
assert(isequaln(MKF_GPB1.x0, zeros(n, 1)))
assert(isequaln(MKF_GPB1.y0, zeros(ny, 1)))
assert(isequal(MKF_GPB1.p_seq_g_Yk_init, [0.5 0.5]'))
assert(isequaln(MKF_GPB1.Skf, nan(n, ny, MKF_GPB1.n_filt)));
assert(isequaln(MKF_GPB1.Kf, nan(n, ny, MKF_GPB1.n_filt)));
assert(isequaln(MKF_GPB1.Pkf, nan(n, n, MKF_GPB1.n_filt)));
assert(isequaln(MKF_GPB1.Xkf_est, nan(n, 1, MKF_GPB1.n_filt)));
assert(isequaln(MKF_GPB1.Ykf_est, nan(ny, 1, MKF_GPB1.n_filt)));
assert(isequal(MKF_GPB1.Pkp1f, repmat(P0, 1, 1, MKF_GPB1.n_filt)));
assert(isequal(MKF_GPB1.Xkp1f_est, zeros(n, 1, MKF_GPB1.n_filt)));
assert(isequal(MKF_GPB1.Ykp1f_est, zeros(ny, 1, MKF_GPB1.n_filt)));

% Redefine this time with initial conditions
MKF_GPB1x0 = MKFObserverGPB1(models,Ts,P0,T,"MKF_GPB1x0",x0);
assert(isequaln(MKF_GPB1x0.x0, x0))
assert(isequal(MKF_GPB1x0.xkp1_est, x0))
assert(isequal(MKF_GPB1x0.Pkp1, P0))
assert(isequal(MKF_GPB1x0.ykp1_est, C1 * x0))
assert(isequal(MKF_GPB1x0.p_seq_g_Yk_init, [0.5; 0.5]))
MKF_GPB1x0y0 = MKFObserverGPB1(models,Ts,P0,T,'MKF_AMM1x0',x0,y0);
assert(isequaln(MKF_GPB1x0y0.x0, x0))
assert(isequaln(MKF_GPB1x0y0.y0, y0))
assert(isequal(MKF_GPB1x0y0.xkp1_est, x0))
assert(isequal(MKF_GPB1x0y0.ykp1_est, y0))

% Also with initial prior probability values
p_seq_g_Yk_init = [0.6; 0.4];
MKF_GPB1x0y0p0 = MKFObserverGPB1(models,Ts,P0,T,"MKF_GPB1x0",x0,y0, ...
    p_seq_g_Yk_init);
assert(isequaln(MKF_GPB1x0y0p0.x0, x0))
assert(isequaln(MKF_GPB1x0y0p0.y0, y0))
assert(isequaln(MKF_GPB1x0y0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))
assert(isequal(MKF_GPB1x0y0p0.xkp1_est, x0))
assert(isequal(MKF_GPB1x0y0p0.ykp1_est, y0))
assert(isequal(MKF_GPB1x0y0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))

% Define GPB2
MKF_GPB2 = MKFObserverGPB2(models,Ts,P0,T,"MKF_GPB2");

assert(strcmp(MKF_GPB2.type, "MKF_GPB2"))
assert(isequal(MKF_GPB2.models, models))
assert(isequal(MKF_GPB2.Ts, Ts))
assert(isequal(MKF_GPB2.P0, P0))
assert(isequal(MKF_GPB2.T, T))
assert(strcmp(MKF_GPB2.label, "MKF_GPB2"))
assert(MKF_GPB2.n_filt == 4)
assert(MKF_GPB2.n == n)
assert(MKF_GPB2.nu == nu)
assert(MKF_GPB2.ny == ny)
assert(MKF_GPB2.nj == 2)
assert(isequal(MKF_GPB2.T, T))
assert(isequal(MKF_GPB2.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_GPB2.Pkp1, P0))
assert(isequal(MKF_GPB2.ykp1_est, 0))
assert(isequaln(MKF_GPB2.xk_est, nan(n, 1)))
assert(isequaln(MKF_GPB2.Pk, nan(n)))
assert(isequaln(MKF_GPB2.yk_est, nan(ny, 1)))
assert(isequaln(MKF_GPB2.gamma_km1, [0 1 0 1]'))
assert(isequaln(MKF_GPB2.gamma_k, [0 0 1 1]'))
assert(isequaln(MKF_GPB2.p_yk_g_seq_Ykm1, nan(MKF_GPB2.n_filt, 1)))
assert(isequaln(MKF_GPB2.p_gammak_g_Ykm1, nan(MKF_GPB2.n_filt, 1)))
assert(isequaln(MKF_GPB2.p_gamma_k_g_gamma_km1, [0.95 0.01 0.05 0.99]'))
assert(isequaln(MKF_GPB2.p_seq_g_Ykm1, nan(MKF_GPB2.n_filt, 1)))
assert(isequaln(MKF_GPB2.x0, zeros(n, 1)))
assert(isequaln(MKF_GPB2.y0, zeros(ny, 1)))
assert(isequal(MKF_GPB2.p_seq_g_Yk_init, ones(4, 1) ./ 4))
assert(isequaln(MKF_GPB2.Skf, nan(n, ny, MKF_GPB2.n_filt)));
assert(isequaln(MKF_GPB2.Kf, nan(n, ny, MKF_GPB2.n_filt)));
assert(isequaln(MKF_GPB2.Pkf, nan(n, n, MKF_GPB2.n_filt)));
assert(isequaln(MKF_GPB2.Xkf_est, nan(n, 1, MKF_GPB2.n_filt)));
assert(isequaln(MKF_GPB2.Ykf_est, nan(ny, 1, MKF_GPB2.n_filt)));
assert(isequal(MKF_GPB2.Pkp1f, repmat(P0, 1, 1, MKF_GPB2.n_filt)));
assert(isequal(MKF_GPB2.Xkp1f_est, zeros(n, 1, MKF_GPB2.n_filt)));
assert(isequal(MKF_GPB2.Ykp1f_est, zeros(ny, 1, MKF_GPB2.n_filt)));

% Choose observers to include in simulation
observers = {KF1, KF2, MKF_AMM1, MKF_GPB1, SKF};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Identify which observer to log MKF data for
f_mkf = 4;

% Simulate observers - without measurement noise (Y)
[Xk_est,Yk_est,DiagPk,Xkp1_est,Ykp1_est,DiagPkp1,MKF_K_obs,MKF_trP_obs, ...
    MKF_i,MKF_p_seq_g_Yk] = run_simulation_obs(Y,U,observers,f_mkf);

% Output estimation errors
E_obs = Y - Yk_est;

% Combine and display results
% sim_results1 = table(t,Gamma,U,X,Y,Ym,Xk_est,Yk_est,E_obs)
% mkf_sim_results = table(t,Gamma,U,X,Y,Ym,MKF_i,MKF_p_seq_g_Yk,MKF_trP_obs)

% figure(2); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,obs_labels)


% Convert to table
E_obs = array2table(Y - Yk_est, 'VariableNames', obs_labels);
Xk_est = array2table(Xk_est, 'VariableNames', obs_labels);

% Compute mean-squared error
rmses = table(sqrt(mean(sum(E_obs{:, :}.^2), 1))', 'RowNames', obs_labels, ...
    'VariableNames', {'RMSE'});

% Check KF1 was accurate before system switched
assert(max(abs(E_obs{t <= 9.5, "KF1"})) < 1e-5)

% Check KF2 was accurate after system switched
assert(max(abs(E_obs{t > 20, "KF2"})) < 1e-3)

% Check SKF matches KF1 before system switched
assert(isequal(Xk_est{t <= 9.5, "KF1"}, Xk_est{t <= 9.5, "SKF1"}))

% Check SKF converges to KF2 after system switched
assert(max(abs(Xk_est{t > 20, "KF2"} - Xk_est{t > 20, "SKF1"})) < 0.01)

% Check MKF_AMM converged to KF1 and remained there
% TODO: Is this result expected for AMM?
assert(max(abs(Xk_est{t > 9.5, "KF1"} - Xk_est{t > 9.5, "MKF_AMM1"})) < 1e-8)

% Check MKF_AMM converged to KF2 after system switched - it does not.
%assert(max(abs(Xk_est{t > 20, "KF2"} - Xk_est{t > 20, "MKF_AMM1"})) < 1e-10)

% Check GPB1
if any(strcmp("MKF_GPB1", Xk_est.Properties.VariableNames))
    % Check GPB1 is close to KF1 before system switched
    assert(abs(Xk_est{t == 9.5, "KF1"} - Xk_est{t == 9.5, "MKF_GPB1"}) < 1e-3)
    % Check GPB1 is close to KF2 after system switched
    assert(max(abs(Xk_est{t > 20, "KF2"} - Xk_est{t > 20, "MKF_GPB1"})) < 0.002)
end

% Check root-mean-squared errors
rmse_test = table([15.2373 4.0535 15.2373 0.2419 0.0000]', ...
    'RowNames', ["KF1"  "KF2" "MKF_AMM1" "MKF_GPB1" "SKF1"], ...
    'VariableNames', {'RMSE'});
for i = 1:size(rmses, 1)
    label = rmses.Properties.RowNames{i};
    assert(round(rmses{label, 'RMSE'}, 4) == rmse_test{label, 'RMSE'})
end

% Test reset methods
% Should set observer variables to original initial conditions
KF1.reset()
KF2.reset()
MKF_AMM1.reset()
MKF_GPB1.reset()
SKF.reset()

assert(isequal(MKF_AMM1.xkp1_est, MKF_AMM1.x0))
assert(isequal(MKF_AMM1.Pkp1, MKF_AMM1.P0))
assert(isequal(MKF_AMM1.ykp1_est, C{1} * MKF_AMM1.x0))
assert(isequal(MKF_AMM1.p_seq_g_Yk, MKF_AMM1.p_seq_g_Yk_init))

assert(isequal(MKF_GPB1.xkp1_est, MKF_GPB1.x0))
assert(isequal(MKF_GPB1.Pkp1, MKF_GPB1.P0))
assert(isequal(MKF_GPB1.ykp1_est, C{1} * MKF_GPB1.x0))
assert(isequal(MKF_GPB1.p_seq_g_Yk, MKF_GPB1.p_seq_g_Yk_init))
assert(isequaln(MKF_GPB1.p_yk_g_seq_Ykm1, nan(MKF_GPB1.n_filt, 1)))
assert(isequaln(MKF_GPB1.p_gammak_g_Ykm1, nan(MKF_GPB1.n_filt, 1)))
assert(isequaln(MKF_GPB1.p_gamma_k_g_gamma_km1, [0.95 0.99]'))  % TODO: Is this correct?
assert(isequaln(MKF_GPB1.p_seq_g_Ykm1, nan(MKF_GPB1.n_filt, 1)))

% Redefine a new observer (identical to above)
MKF_GPB1_new = MKFObserverGPB1(models,Ts,P0,T,'MKF_GPB1');
assert(isequaln(MKF_GPB1_new, MKF_GPB1))
MKF_GPB1_new.label = "MKF_GPB1_new";

% Make a copy
MKF_GPB1_copy = MKF_GPB1_new.copy();
assert(isequaln(MKF_GPB1_copy, MKF_GPB1_new))
MKF_GPB1_copy.label = "MKF_GPB1_copy";

% Choose observers to include in simulation
observers = {KF1, KF2, SKF, MKF_GPB1, MKF_GPB1_new, MKF_GPB1_copy};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Identify which observer to log MKF data for
f_mkf = 4;

% Simulate observers - with measurement noise (Ym)
[Xk_est,Yk_est,DiagPk,Xkp1_est,Ykp1_est,DiagPkp1,MKF_K_obs,MKF_trP_obs, ...
    MKF_i,MKF_p_seq_g_Yk] = run_simulation_obs(Y,U,observers,f_mkf);

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

% Check final state estimates
test_X_est = [-1.191082  9.901224  9.901227  9.901228  9.901228  9.901228];
assert(isequal(round(Xk_est(t == t(end), :), 6), test_X_est))

% Check final error covariance estimates
test_DiagP = [0.015011  0.022516  0.022516  0.022516  0.022516  0.022516];
assert(isequal(round(DiagPk(t == t(end), :), 6), test_DiagP))

% Compute mean-squared error
mses = nanmean(E_obs.^2);
%array2table(mses,'VariableNames',obs_labels)

% Check MKF observer estimation error
assert(round(mses(f_mkf), 6) == 0.000959)

% Check all observer estimation errors
assert(isequal(round(mses, 6), ...
    [3.806151 0.269363 0 0.000959 0.000959 0.000959]))

% % Plot selected observers
% figure(4); clf
% plot_obs_estimates(t,X,X_est(:,[3 4]),Y,Y_est(:,[3 4]),obs_labels([3 4]))


%% Simulation test on 2x2 system with random shock

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
C = [ 0.1110       0  0  0;
           0  0.1110  0  0];
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
nT = 100;
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
Gamma = zeros(nT+1, 2);
Gamma(t == t_shock(1), 1) = 1;
Gamma(t == t_shock(2), 2) = 1;
Wp = du0' .* Gamma;

U_sim = [U Wp];

% Observer parameters (same for all observers)
P0 = 1000*eye(n);
x0 = [0.5 -0.1 0.1 -0.1]';  % optional
y0 = C * x0;  % optional

% 4 possible process noise covariance matrices
Qj = {diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,1)^2]), ...
      diag([0.01 0.01 sigma_wp(1,2)^2 sigma_wp(2,1)^2]), ...
      diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,2)^2]), ...
      diag([0.01 0.01 sigma_wp(1,2)^2 sigma_wp(2,2)^2])};
R = diag(sigma_M.^2);

% 3 system models
nj = 3;
m1.A = A;
m1.B = Bu;
m1.C = C;
m1.Q = Qj{1};
m1.R = R;
m2.A = A;
m2.B = Bu;
m2.C = C;
m2.Q = Qj{2};
m2.R = R;
m3.A = A;
m3.B = Bu;
m3.C = C;
m3.Q = Qj{3};
m3.R = R;
models = {m1, m2, m3};

% Custom MKF observer with correct shock sequence
seq = {zeros(1, nT+1); zeros(1, nT+1); zeros(1, nT+1); zeros(1, nT+1)};
seq{2}(t == t_shock(1)) = 1;  % shock 1
seq{3}(t == t_shock(2)) = 2;  % shock 2
seq{4}(t == t_shock(1)) = 1;  % both
seq{4}(t == t_shock(2)) = 2;
p_gamma = [1-epsilon epsilon]';
Z = [0 0; 1 0; 0 1];  % combinations
p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
p_gamma = p_gamma ./ sum(p_gamma);  % normalized
T = repmat(p_gamma', nj, 1);

% Define MKF observers

% 1. MKF_AMM
MKF_AMM1 = MKFObserverAMM(models,Ts,P0,'MKF_AMM1');

assert(strcmp(MKF_AMM1.type, "MKF_AMM"))
assert(isequal(MKF_AMM1.models, models))
assert(isequal(MKF_AMM1.Ts, Ts))
assert(isequal(MKF_AMM1.P0, P0))
assert(strcmp(MKF_AMM1.label, "MKF_AMM1"))
assert(MKF_AMM1.n_filt == 3)
assert(MKF_AMM1.n == n)
assert(MKF_AMM1.nu == nu)
assert(MKF_AMM1.ny == ny)
assert(MKF_AMM1.nj == 3)
assert(isequal(MKF_AMM1.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_AMM1.Pkp1, P0))
assert(isequal(MKF_AMM1.ykp1_est, zeros(ny, 1)))
assert(isequaln(MKF_AMM1.xk_est, nan(n, 1)))
assert(isequaln(MKF_AMM1.Pk, nan(n)))
assert(isequaln(MKF_AMM1.yk_est, nan(ny, 1)))
assert(isequaln(MKF_AMM1.x0, zeros(n, 1)))
assert(isequaln(MKF_AMM1.y0, zeros(ny, 1)))
assert(isequal(MKF_AMM1.p_seq_g_Yk_init, ones(3, 1) ./ 3))

% Redefine this time with initial conditions
MKF_AMM1x0 = MKFObserverAMM(models,Ts,P0,'MKF_AMM1x0',x0);
assert(isequaln(MKF_AMM1x0.x0, x0))
assert(isequal(MKF_AMM1x0.xkp1_est, x0))
assert(isequal(MKF_AMM1x0.Pkp1, P0))
assert(isequal(MKF_AMM1x0.ykp1_est, C * x0))
assert(isequal(MKF_AMM1x0.p_seq_g_Yk_init, ones(3, 1) ./ 3))
MKF_AMM1x0y0 = MKFObserverAMM(models,Ts,P0,'MKF_AMM1x0',x0,y0);
assert(isequaln(MKF_AMM1x0y0.x0, x0))
assert(isequaln(MKF_AMM1x0y0.y0, y0))
assert(isequal(MKF_AMM1x0y0.xkp1_est, x0))
assert(isequal(MKF_AMM1x0y0.ykp1_est, y0))

% Also with initial prior probability values
p_seq_g_Yk_init = [0.6; 0.4];
MKF_AMM1x0y0p0 = MKFObserverAMM(models,Ts,P0,"MKF_GPB1x0",x0,y0, ...
    p_seq_g_Yk_init);
assert(isequaln(MKF_AMM1x0y0p0.x0, x0))
assert(isequaln(MKF_AMM1x0y0p0.y0, y0))
assert(isequaln(MKF_AMM1x0y0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))
assert(isequal(MKF_AMM1x0y0p0.xkp1_est, x0))
assert(isequal(MKF_AMM1x0y0p0.ykp1_est, y0))
assert(isequal(MKF_AMM1x0y0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))

% 2. MKF_GPB1
MKF_GPB1 = MKFObserverGPB1(models,Ts,P0,T,'MKF_GPB1');

assert(strcmp(MKF_GPB1.type, "MKF_GPB1"))
assert(isequal(MKF_GPB1.models, models))
assert(isequal(MKF_GPB1.Ts, Ts))
assert(isequal(MKF_GPB1.P0, P0))
assert(isequal(MKF_GPB1.T, T))
assert(strcmp(MKF_GPB1.label, "MKF_GPB1"))
assert(MKF_GPB1.n_filt == 3)
assert(MKF_GPB1.n == n)
assert(MKF_GPB1.nu == nu)
assert(MKF_GPB1.ny == ny)
assert(MKF_GPB1.nj == 3)
assert(isequal(MKF_GPB1.T, T))
assert(isequal(MKF_GPB1.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_GPB1.Pkp1, P0))
assert(isequal(MKF_GPB1.ykp1_est, zeros(ny, 1)))
assert(isequaln(MKF_GPB1.xk_est, nan(n, 1)))
assert(isequaln(MKF_GPB1.Pk, nan(n)))
assert(isequaln(MKF_GPB1.yk_est, nan(ny, 1)))
assert(isequaln(MKF_GPB1.p_yk_g_seq_Ykm1, nan(MKF_GPB1.n_filt, 1)))
assert(isequaln(MKF_GPB1.p_gammak_g_Ykm1, nan(MKF_GPB1.n_filt, 1)))
assert(isequaln(round(MKF_GPB1.p_gamma_k_g_gamma_km1, 6), ...
    [0.980198 0.009901 0.009901]'))  % Doesn this make sense?
assert(isequaln(MKF_GPB1.p_seq_g_Ykm1, nan(MKF_GPB1.n_filt, 1)))
assert(isequaln(MKF_GPB1.x0, zeros(n, 1)))
assert(isequaln(MKF_GPB1.y0, zeros(ny, 1)))
assert(isequal(MKF_GPB1.p_seq_g_Yk_init, ones(3, 1) ./ 3))
assert(isequaln(MKF_GPB1.Skf, nan(n, ny, MKF_GPB1.n_filt)));
assert(isequaln(MKF_GPB1.Kf, nan(n, ny, MKF_GPB1.n_filt)));
assert(isequaln(MKF_GPB1.Pkf, nan(n, n, MKF_GPB1.n_filt)));
assert(isequaln(MKF_GPB1.Xkf_est, nan(n, 1, MKF_GPB1.n_filt)));
assert(isequaln(MKF_GPB1.Ykf_est, nan(ny, 1, MKF_GPB1.n_filt)));
assert(isequal(MKF_GPB1.Pkp1f, repmat(P0, 1, 1, MKF_GPB1.n_filt)));
assert(isequal(MKF_GPB1.Xkp1f_est, zeros(n, 1, MKF_GPB1.n_filt)));
assert(isequal(MKF_GPB1.Ykp1f_est, zeros(ny, 1, MKF_GPB1.n_filt)));

% Redefine this time with initial conditions
MKF_GPB1x0 = MKFObserverGPB1(models,Ts,P0,T,"MKF_GPB1x0",x0);
assert(isequaln(MKF_GPB1x0.x0, x0))
assert(isequal(MKF_GPB1x0.xkp1_est, x0))
assert(isequal(MKF_GPB1x0.Pkp1, P0))
assert(isequal(MKF_GPB1x0.ykp1_est, C * x0))
assert(isequal(MKF_GPB1x0.p_seq_g_Yk_init, ones(3, 1) ./ 3))
MKF_GPB1x0y0 = MKFObserverGPB1(models,Ts,P0,T,'MKF_AMM1x0',x0,y0);
assert(isequaln(MKF_GPB1x0y0.x0, x0))
assert(isequaln(MKF_GPB1x0y0.y0, y0))
assert(isequal(MKF_GPB1x0y0.xkp1_est, x0))
assert(isequal(MKF_GPB1x0y0.ykp1_est, y0))

% Also with initial prior probability values
p_seq_g_Yk_init = [0.6; 0.4];
MKF_GPB1x0y0p0 = MKFObserverGPB1(models,Ts,P0,T,"MKF_GPB1x0",x0,y0, ...
    p_seq_g_Yk_init);
assert(isequaln(MKF_GPB1x0y0p0.x0, x0))
assert(isequaln(MKF_GPB1x0y0p0.y0, y0))
assert(isequaln(MKF_GPB1x0y0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))
assert(isequal(MKF_GPB1x0y0p0.xkp1_est, x0))
assert(isequal(MKF_GPB1x0y0p0.ykp1_est, y0))
assert(isequal(MKF_GPB1x0y0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))

% Define scheduled Kalman filter
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
% seq = sum(alpha .* 2.^(1:-1:0), 2)';
Aj = repmat({A}, 1, nj);
Buj = repmat({Bu}, 1, nj);
Cj = repmat({C}, 1, nj);
Rj = repmat({R}, 1, nj);
SKF = MKFObserverSchedF(Aj,Buj,Cj,Ts,P0,Qj,Rj,seq{4},"SKF");

% Choose observers to test
observers = {MKF_AMM1, MKF_GPB1, SKF};
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
Ym = Y + sigma_MP'.*randn(nT+1, ny);

% Identify which observer to log MKF data for
f_mkf = 1;

% Simulate observers
[Xk_est,Yk_est,DiagPk,Xkp1_est,Ykp1_est,DiagPkp1,MKF_K_obs,MKF_trP_obs, ...
    MKF_i,MKF_p_seq_g_Yk] = run_simulation_obs(Ym,U,observers,f_mkf);

% Output estimation errors
E_obs = repmat(Y, 1, n_obs) - Yk_est;

% Combine and display results
sim_results1 = table(t,Gamma,U,X,Y,Ym,Xk_est,Yk_est,E_obs);

% Plot observer estimates
% figure(5); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,obs_labels)

% Average results and convert to table
E_obs_avg = array2table(reshape(mean(reshape(E_obs,[],ny,n_obs), 2),[],n_obs), ...
    'VariableNames', obs_labels);
Xk_est_avg = array2table(reshape(mean(reshape(Xk_est,[],n,n_obs), 2),[],n_obs), ...
    'VariableNames', obs_labels);
DiagPk_avg = array2table(reshape(mean(reshape(DiagPk,[],n,n_obs), 2),[],n_obs), ...
    'VariableNames', obs_labels);

% Check final state estimates
test_X_est_avg = [2.293237    2.301738    2.301600];
assert(isequal(round(Xk_est_avg{t == t(end), :}, 6), test_X_est_avg))

% Check final error covariance estimates
% TODO: Haven't checked if these are correct.
test_DiagPk_avg = [0.0427    0.1547    0.0427];
assert(isequal(round(DiagPk_avg{t == t(end), :}, 4), test_DiagPk_avg))

% Compute mean-squared error
rmses = array2table(reshape(sqrt(mean(reshape(E_obs,[],ny,n_obs).^2, 2)),[],n_obs), ...
    'VariableNames', obs_labels);
rmses = array2table(round(mean(rmses{:,:}, 1), 6), 'VariableNames', obs_labels);
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

% Results
MSE_test_values = array2table(...
    [0.025275 0.010339 0.002565], ...
    'VariableNames', {'MKF_AMM1', 'MKF_GPB1', 'SKF'});
assert(isequal(rmses, MSE_test_values))



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
%Gpss1 = ss(A1,B1,C1,D1,Ts);

% Model #2
A2 = 0.9;
B2 = 1;
C2 = -0.3;  % -ve gain!
D2 = 0;
%Gpss2 = ss(A2,B2,C2,D2,Ts);

% Dimensions
n = size(A1, 1);
nu = size(B1, 2);
ny = size(C1, 1);

% Check dimensions
assert(isequal(size(A1), size(A2)))
assert(isequal(size(B1), size(B2)))
assert(isequal(size(C1), size(C2)))
assert(isequal(size(D1), size(D2)))

% % Define system models
% A = {A1, A2};
% B = {B1, B2};
% C = {C1, C2};
% D = {D1, D2};

% Observer parameters (same for all observers)
P0 = 10000;
x0 = 0.5;  % optional
y0 = 0.1;  % optional
Q1 = 0.01;
R1 = 0.1^2;
Q2 = 0.01;
R2 = 0.1^2;

% Switching parameters
% Q = {Q1,Q2};
% R = {R1,R2};

% Transition probabilities
epsilon = 0.05;
T = [1-epsilon epsilon; epsilon 1-epsilon];
assert(all(sum(T, 2) == 1))

% System indicator sequences
nT = 100;
seq1 = {
    zeros(1, nT+1);
    [zeros(1, 20) ones(1, nT+1-20)];  % equal to Gamma'
    [zeros(1, 40) ones(1, nT+1-40)];
    ones(1, nT+1);
 };

% Define MKF observers

% Two system models
m1.A = A1;
m1.B = B1;
m1.C = C1;
m1.Q = Q1;
m1.R = R1;
m2.A = A2;
m2.B = B2;
m2.C = C2;
m2.Q = Q2;
m2.R = R2;
models = {m1, m2};

% Define multi-model observer with initial conditions
MKF = MKFObserverGPB1(models,Ts,P0,T,'MKF',x0);

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

MKF.label = "New name";
assert(~isequal(MKF_copy.label, "New name"))

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


function [Xk_est,Yk_est,DiagPk,Xkp1_est,Ykp1_est,DiagPkp1,MKF_K_obs,MKF_trP_obs, ...
    MKF_i,MKF_p_seq_g_Yk] = run_simulation_obs(Ym,U,observers,f_mkf)
% Simulate observers

    nT = size(Ym, 1) - 1;
    ny = size(Ym, 2);
    n_obs = numel(observers);
    n = size(observers{1}.xkp1_est, 1);

    obs_mkf = observers{f_mkf};
    n_filters = obs_mkf.n_filt;

    Xk_est = zeros(nT+1, n*n_obs);
    Yk_est = zeros(nT+1, ny*n_obs);
    Xkp1_est = zeros(nT+1, n*n_obs);
    Ykp1_est = zeros(nT+1, ny*n_obs);
    DiagPk = nan(nT+1, n*n_obs);
    DiagPkp1 = nan(nT+1, n*n_obs);
    MKF_K_obs = cell(nT+1, n*n_filters);
    MKF_trP_obs = nan(nT+1, n_filters);
    MKF_i = nan(nT+1, 2);
    MKF_p_seq_g_Yk = nan(nT+1, n_filters);

    for i = 1:nT+1

        yk = Ym(i, :)';
        uk = U(i, :)';

        % Arrays to store results
        xk_est = nan(1, n*n_obs);
        yk_est = nan(1, ny*n_obs);
        xkp1_est = nan(1, n*n_obs);
        ykp1_est = nan(1, ny*n_obs);
        diagPk = nan(1, n*n_obs);
        diagPkp1 = nan(1, n*n_obs);

        % Update observers
        for f = 1:n_obs
            obs = observers{f};
            obs.update(yk, uk);
            if f == f_mkf
                if isprop(obs, "i")
                    MKF_i(i, :) = obs.i;
                else
                    MKF_i(i, :) = nan;
                end
                MKF_p_seq_g_Yk(i, :) = obs.p_seq_g_Yk';
                for j = 1:obs.n_filt
                    switch obs.type
                        case {"KF", "SKF", "MKF"}
                            % Only works for ny = 1
                            %MKF_K_obs{i, j} = obs.filters{j}.K';
                            MKF_trP_obs(i, j) = trace(obs.filters{j}.P);
                        case {"KFF", "SKFF"}
                            %MKF_K_obs{i, j} = obs.filters{j}.Kf';
                            MKF_trP_obs(i, j) = trace(obs.filters{j}.Pkp1);
                        case {"MKF_AMM", "MKF_GPB1"}
                            MKF_trP_obs(i, j) = trace(obs.Pkf(:, :, j));
                    end
                end
                switch obs.type
                end
            end
            if isprop(obs,'xk_est')
                xk_est(1, (f-1)*n+1:f*n) = obs.xk_est';
                yk_est(1, (f-1)*ny+1:f*ny) = obs.yk_est';
                diagPk(1, (f-1)*n+1:f*n) = diag(obs.Pk)';
            end
            if isprop(obs,'xkp1_est')
                xkp1_est(1, (f-1)*n+1:f*n) = obs.xkp1_est';
                ykp1_est(1, (f-1)*ny+1:f*ny) = obs.ykp1_est';
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
