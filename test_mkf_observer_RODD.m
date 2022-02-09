% Test functions mkf_observer_RODD.m and update_MKF.m

clear all
plot_dir = 'plots';

seed = 0;
rng(seed)


%% Test observers for SISO system

% Load system and disturbance model from file
sys_rodin_step

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);
nw = sum(~u_meas);

% Set noise variances for observer design
sigma_M = 0.1;
sigma_W = [0; 0];

% Load observers from file
obs_rodin_step

% Check observer attributes
assert(MKF1.epsilon == epsilon)
assert(isequal(MKF1.sigma_wp, sigma_wp))
assert(MKF1.n_filt == 7)
assert(isequaln(MKF1.i, [0 0]))
assert(MKF1.n == 2)
assert(MKF1.nu == 1)
assert(MKF1.ny == 1)
assert(MKF1.nj == 2)
assert(isequal(MKF1.A{1}, A) && isequal(MKF1.A{2}, A))
assert(isequal(MKF1.B{1}, Bu) && isequal(MKF1.B{2}, Bu))
assert(isequal(MKF1.C{1}, C) && isequal(MKF1.C{2}, C))
assert(isequal(MKF1.D{1}, Du) && isequal(MKF1.D{2}, Du))
assert(MKF1.Ts == Ts)
assert(isequaln(MKF1.u_meas, u_meas))
assert(isequal(size(MKF1.Q), [1 2]))
assert(isequal(MKF1.Q{1}, [0.01 0; 0 sigma_wp(1)^2/MKF1.d]))
assert(isequal(MKF1.Q{2}, [0.01 0; 0 sigma_wp(2)^2/MKF1.d]))
assert(isequal(size(MKF1.R), [1 2]))
assert(isequal(MKF1.R{1}, R) && isequal(MKF1.R{2}, R))
assert(numel(MKF1.filters) == MKF1.n_filt)
assert(isequal(size(MKF1.seq), [MKF1.n_filt 1]))
assert(isequal(size(cell2mat(MKF1.seq)), [MKF1.n_filt MKF1.f]))
assert(MKF1.beta == sum(MKF1.p_seq))
assert(MKF1.f == size(MKF1.seq{1}, 2))
assert(isequal(MKF1.xkp1_est, zeros(n,1)))
assert(isequal(MKF1.ykp1_est, zeros(ny,1)))
assert(isequal(round(MKF1.alpha, 4), 0.0490))
assert(isequal(round(MKF1.p_gamma, 4), [0.9510; 0.0490]))

assert(MKF2.epsilon == epsilon)
assert(isequal(MKF2.sigma_wp, sigma_wp))
assert(MKF2.n_filt == 56)
assert(isequaln(MKF1.i, [0 0]))
assert(MKF2.n == 2)
assert(MKF2.nu == 1)
assert(MKF2.ny == 1)
assert(MKF2.nj == 2)
assert(isequal(MKF2.A{1}, A) && isequal(MKF2.A{2}, A))
assert(isequal(MKF2.B{1}, Bu) && isequal(MKF2.B{2}, Bu))
assert(isequal(MKF2.C{1}, C) && isequal(MKF2.C{2}, C))
assert(isequal(MKF2.D{1}, Du) && isequal(MKF2.D{2}, Du))
assert(MKF2.Ts == Ts)
assert(isequaln(MKF2.u_meas, u_meas))
assert(isequal(size(MKF2.Q), [1 2]))
assert(isequal(MKF2.Q{1}, [0.01 0; 0 sigma_wp(1)^2/MKF2.d]))
assert(isequal(MKF2.Q{2}, [0.01 0; 0 sigma_wp(2)^2/MKF2.d]))
assert(isequal(size(MKF1.R), [1 2]))
assert(isequal(MKF2.R{1}, R) && isequal(MKF2.R{2}, R))
assert(numel(MKF2.filters) == MKF2.n_filt)
assert(isequal(size(MKF2.seq), [MKF2.n_filt 1]))
assert(isequal(size(cell2mat(MKF2.seq)), [MKF2.n_filt MKF2.f]))
assert(MKF2.beta == sum(MKF2.p_seq))
assert(MKF2.f == size(MKF2.seq{1}, 2))
assert(isequal(MKF2.xkp1_est, zeros(n,1)))
assert(isequal(MKF2.ykp1_est, zeros(ny,1)))
assert(isequal(round(MKF2.alpha, 4), MKF2.epsilon))
assert(isequal(round(MKF2.p_gamma, 4), [1-MKF2.epsilon; MKF2.epsilon]))

% Check optional definition with an initial state estimate works
x0 = [0.1; 0.5];
MKF_testx0 = mkf_observer_RODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label,x0);
assert(isequal(MKF_testx0.xkp1_est, x0))
assert(isequal(MKF_testx0.ykp1_est, C * x0))

% Simulation settings
nT = 100;
t = Ts*(0:nT)';

% Choose time and amplitude of input disturbance
t_shock = 9.5;
du0 = 1;
% When you make the shock larger the MKF observers
% do better
%du0 = 2;

% Measured input
%U = (idinput(size(t)) + 1)/2;
U = zeros(size(t));
U(t >= 1) = -1;

% Disturbance input
alpha = zeros(nT+1, 1);
alpha(t == t_shock) = 1;  % this is used by the SKF observer
Wp = du0 .* alpha;

U_sim = [U Wp];

% Custom MKF test observer

% Devise a custom multi-model filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% this test simulation (t == t_shock)

% Multiple model filter - two sequences, one empty, one correct
A2 = repmat({A}, 1, 2);
Bu2 = repmat({Bu}, 1, 2);
C2 = repmat({C}, 1, 2);
Du2 = repmat({Du}, 1, 2);
P0 = 1000*eye(n);
Q0 = diag([0.01 1]);
P0_init = repmat({P0}, 1, 2);
Q2 = {diag([Q0(1,1) sigma_wp(1,1)^2]), ...
      diag([Q0(1,1) sigma_wp(1,2)^2])};
R2 = {sigma_M.^2, sigma_M.^2};
seq = {zeros(1, nT+1); zeros(1, nT+1)};
seq{2}(t == t_shock) = 1;
p_gamma = [1-epsilon epsilon]';
T = repmat(p_gamma', 2, 1);
d = 1;
MKF3 = mkf_observer(A2,Bu2,C2,Du2,Ts,P0_init,Q2,R2,seq,T,d,'MKF3');

% Multiple model filter - one sequence with correct shock
seq = {zeros(1, nT+1)};
seq{1}(t == t_shock) = 1;
p_gamma = [1-epsilon epsilon]';
T = repmat(p_gamma', 2, 1);
d = 1;
MKF4 = mkf_observer(A2,Bu2,C2,Du2,Ts,P0_init,Q2,R2,seq,T,d,'MKF4');


% Choose observers to test
observers = {KF2, KF3, SKF, MKF1, MKF2, MKF3, MKF4};


% Note: KF1 is too slow to pass static error test here

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
sigma_MP = 0;  % Set to zero for testing
Y_m = Y + sigma_MP'.*randn(size(Y));


% Simulate observers

% Measured inputs (not including disturbances)
U_m = U;

n_obs = numel(observers);
MSE = containers.Map();
for i = 1:n_obs

    obs = observers{i};
    [obs, sim_results] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m,obs, ...
        alpha,u_meas);

    % Check observer errors are zero prior to
    % input disturbance
    assert(all(abs(sim_results.X_est(1:20,:) - X(1:20, :)) < 1e-10, [1 2]))
    assert(all(abs(sim_results.Y_est(1:20,:) - Y(1:20, :)) < 1e-10))

    % Check observer static errors are small
    % after input disturbance
    % TODO: Should these be closer?
    if all(sigma_MP == 0)
        assert(abs(sim_results.Y_est(end, :) - Y(end, :)) < 1e-3);
        assert(abs(sim_results.X_est(end, 2) - du0) < 1e-3);
    end

    % Compute mean-squared error
    Y_est = sim_results.Y_est;
    MSE(obs.label) = mean((Y_est - Y).^2);
    %fprintf("%d, %s: %f\n", i, obs.label, mean((Y_est - Y).^2))

    % Save updated observer
    observers{i} = obs;

end


% Display results of last simulation

% X_est = sim_results.X_est;
% E_obs = sim_results.E_obs;
% K_obs = sim_results.K_obs;
% trP_obs = sim_results.trP_obs;
% 
% table(t,alpha,U,Wp,X,Y,Y_m,X_est,Y_est,E_obs)

% % Display gains and trace of covariance matrix
% table(t, cell2mat(K_obs), cell2mat(trP_obs), ...
%     'VariableNames',{'t', 'K{1}, K{2}', 'trace(P{1}), trace(P{2})'})

% % Show table of mean-squared errors
% table(MSE.keys', cell2mat(MSE.values'), ...
%     'VariableNames', {'Observer', 'MSE'})


% % Plot of inputs and outputs
% 
% obs_label = obs.label;
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% figure(1); clf
% colors = get(gca,'colororder');
% ax1 = subplot(4,1,1);
% stairs(t,Y_m); hold on
% stairs(t,Y_est,'Linewidth',2);
% ax1.ColorOrder = colors(1:size(Y_m,2),:);
% max_min = [min(min([Y_m Y_est])) max(max([Y_m Y_est]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$y_m(k)$ and $\hat{y}(k)$')
% title(sprintf('%s - Process output measurements and estimates', obs_label))
% legend('$y_m(k)$','$\hat{y}(k)$')
% grid on
% 
% ax2 = subplot(4,1,2);
% stairs(t,X); hold on
% stairs(t,X_est,'Linewidth',2);
% ax2.ColorOrder = colors(size(Y,2)+1:size(Y,2)+size(X,2),:);
% max_min = [min(min([X X_est])) max(max([X X_est]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$x_i(k)$ and $\hat{x}_i(k)$')
% labels = repmat({''}, 1, n*2);
% for i=1:n
%     labels{i} = sprintf("$x_{%d}(k)$", i);
% end
% for i=1:n
%     labels{i+n} = sprintf("$%s{x}_{%d}(k)$", '\hat', i);
% end
% legend(labels)
% title('States and state estimates')
% grid on
% 
% ax3 = subplot(4,1,3);
% stairs(t,U,'Linewidth',2); hold on;
% stairs(t,Wp,'Linewidth',2)
% max_min = [min(min([U Wp])) max(max([U Wp]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$u(k)$ and $w_p(k)$')
% legend('$u(k)$', '$w_p(k)$')
% title('Process inputs')
% grid on
% 
% ax4 = subplot(4,1,4);
% stairs(t,alpha,'Color',colors(end,:),'Linewidth',2)
% max_min = [min(min(alpha)) max(max(alpha))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$\gamma(k)$')
% title('Random shock sequence')
% grid on
% 
% linkaxes([ax1, ax2 ax3 ax4], 'x')
% 
% set(gcf,'Position',[100 200 560 600]);
% 
% 
% % Plot of conditional filter probabilities
% 
% switch obs.label
%     case {'MKF1', 'MKF2'}
%         p_seq_g_Yk = sim_results.MKF_p_seq_g_Yk;
%         % Note: first data points are nans,
%         % ignore last data point to make plot wider
% 
%         figure(11); clf
%         t = Ts*(0:nT)';
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$Pr(\Gamma(k) \mid Y(k))$'};
%         filename = sprintf('rod_mkf_observer_test_pyk_wfplot.png');
%         filepath = fullfile(plot_dir, filename);
%         show_waterfall_plot(t(2:end-1), p_seq_g_Yk(2:end-1, :), [0 1], ...
%             ax_labels, [0 82], filepath);
%         title('Conditional probabilities of y(k)')
% end
% 
% 
% % Plot of trace of filter covariance matrices
% 
% switch obs.label
%     case {'MKF1', 'MKF2'}
%         trP_obs = cell2mat(sim_results.trP_obs);
% 
%         figure(12); clf
%         t = Ts*(0:nT)';
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$Tr(P(k))$'};
%         filename = sprintf('rod_mkf_observer_test_trP_wfplot.png');
%         filepath = fullfile(plot_dir, filename);
%         show_waterfall_plot(t, trP_obs, [0 5], ax_labels, [0 82], filepath);
%         title('Trace of covariance matrices')
% 
% end
% 
% % Plot of filter correction gains (k1)
% 
% switch obs.label
%     case {'MKF1', 'MKF2'}
%         K_obs = cell2mat(sim_results.K_obs);
%         % Select first gain value onlu
%         K1_obs = K_obs(:,1:2:end);
% 
%         figure(13); clf
%         t = Ts*(0:nT)';
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$Tr(P(k))$'};
%         filename = sprintf('rod_mkf_observer_test_K_wfplot.png');
%         filepath = fullfile(plot_dir, filename);
%         show_waterfall_plot(t, K1_obs, [0 5], ax_labels, [0 82], filepath);
%         title('Filter correction gains (k1)')
% 
% end

% Earlier test results (with a shock of amplitude 1)
% MSE_test_values = containers.Map(...
%     {'KF2',   'KF3',   'MKF1',  'MKF2',  'MKF3',  'MKF4',  'SKF'}, ...
%     [0.000934 0.003524 0.004914 0.005016 0.002709 0.000929 0.000929] ...
% );

% Results on Nov 8 before reverting back the Bayesian updating
%MSE_test_values = containers.Map(...
%  {'KF2',   'KF3',   'MKF1',  'MKF2',  'MKF3',  'MKF4',  'SKF'}, ...
%  [0.000934 0.003524 0.009456 0.005016 0.002709 0.000929 0.000929] ...
%);
% Changes since previous results: 
%  - f, m, d parameters for MK1 changed.
%  - Shock probability and variance modified to reflect detection
%    intervals.
%  - Bayesian prob updates only at end of detection intervals
% Note: performance of MKF3 increases if shock amplitude is increased.

% Results on Nov 8 after reverting back the Bayesian updating
MSE_test_values = containers.Map(...
 {'KF2',   'KF3',   'MKF1',  'MKF2',  'MKF3',  'MKF4',  'SKF'}, ...
 [0.000934 0.003524 0.010423 0.005016 0.002709 0.000929 0.000929] ...
);

for label = MSE.keys
   assert(isequal(round(MSE(label{1}), 6), MSE_test_values(label{1})))
end


%% Test initialization of MKF observers on 2x2 system

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

% Kalman filter 1 - tuned to sigma_wp(1)
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,1)^2]);
R = diag(sigma_M.^2);
KF1 = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,'KF1');

% Kalman filter 2 - tuned to sigma_wp(2)
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([0.01 0.01 sigma_wp(1,2)^2 sigma_wp(2,2)^2]);
R = diag(sigma_M.^2);
KF2 = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,'KF2');

% Kalman filter 3 - manually tuned
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([0.01 0.01 0.1^2 0.1^2]);
R = diag(sigma_M.^2);
KF3 = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,'KF3');

% Scheduled Kalman filter
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
SKF = kalman_filter(A,Bu,C,Du,Ts,P0,Q0,R,'SKF');
SKF.Q0 = Q0;
SKF.Bw = Bw;
SKF.sigma_wp = sigma_wp;

% Multiple model filter 1
label = 'MKF1';
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
f = 3;  % 5 fusion horizon
m = 1;  % 1 maximum number of shocks
d = 2;  % 10 spacing parameter
MKF1 = mkf_observer_RODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Multiple model filter 2
label = 'MKF2';
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
f = 5;  % 10 fusion horizon
m = 2;  % 2 maximum number of shocks
d = 2;  % 5 spacing parameter
MKF2 = mkf_observer_RODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Tests for MKF1
assert(isequal(MKF1.epsilon, epsilon))
assert(isequal(MKF1.sigma_wp, sigma_wp))
assert(MKF1.n_filt == 7)
assert(isequaln(MKF1.i, [0 0]))
assert(MKF1.n == 4)
assert(MKF1.nu == 2)
assert(MKF1.ny == 2)
assert(MKF1.nj == 3)
assert(isequal(MKF1.A{1}, A) && isequal(MKF1.A{2}, A))
assert(isequal(MKF1.B{1}, Bu) && isequal(MKF1.B{2}, Bu))
assert(isequal(MKF1.C{1}, C) && isequal(MKF1.C{2}, C))
assert(isequal(MKF1.D{1}, Du) && isequal(MKF1.D{2}, Du))
assert(MKF1.Ts == Ts)
assert(isequaln(MKF1.u_meas, u_meas))
assert(isequal(size(MKF1.Q), [1 3]))
assert(isequal(MKF1.Q{1}, ...
    diag([0.01 0.01 sigma_wp(1, 1)^2/MKF1.d sigma_wp(2, 1)^2/MKF1.d])))
assert(isequal(MKF1.Q{2}, ...
    diag([0.01 0.01 sigma_wp(1, 1)^2/MKF1.d sigma_wp(2, 2)^2/MKF1.d])))
assert(isequal(MKF1.Q{3}, ...
    diag([0.01 0.01 sigma_wp(1, 2)^2/MKF1.d sigma_wp(2, 1)^2/MKF1.d])))
assert(isequal(size(MKF1.R), [1 3]))
assert(isequal(MKF1.R{1}, R) && isequal(MKF1.R{2}, R) && isequal(MKF1.R{3}, R))
assert(numel(MKF1.filters) == MKF1.n_filt)
assert(isequal(size(MKF1.seq), [MKF1.n_filt 1]))
assert(isequal(size(cell2mat(MKF1.seq)), [MKF1.n_filt MKF1.f]))
assert(MKF1.beta == sum(MKF1.p_seq))
assert(MKF1.f == size(MKF1.seq{1}, 2))
assert(isequal(size(MKF1.xkp1_est), [n 1]))
assert(isequal(size(MKF1.ykp1_est), [ny 1]))
assert(sum(MKF1.p_gamma) == 1)
alpha = (1 - (1 - epsilon).^d);  % prob. over detection interval 
p_gamma = [1-alpha'; alpha'];
Z = [0 0; 0 1; 1 0];  % combinations
p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
p_gamma = p_gamma ./ sum(p_gamma);  % normalized
assert(isequal(round(p_gamma, 6), [0.960977; 0.019512; 0.019512]))
assert(isequal(round(MKF1.p_gamma, 6), [0.960977; 0.019512; 0.019512]))

% Tests for MKF2
assert(isequal(MKF2.epsilon, epsilon))
assert(isequal(MKF2.sigma_wp, sigma_wp))
assert(MKF2.n_filt == 56)
assert(isequaln(MKF2.i, [0 0]))
assert(MKF2.n == 4)
assert(MKF2.nu == 2)
assert(MKF2.ny == 2)
assert(MKF2.nj == 4)
assert(isequal(MKF2.A{1}, A) && isequal(MKF2.A{2}, A))
assert(isequal(MKF2.B{1}, Bu) && isequal(MKF2.B{2}, Bu))
assert(isequal(MKF2.C{1}, C) && isequal(MKF2.C{2}, C))
assert(isequal(MKF2.D{1}, Du) && isequal(MKF2.D{2}, Du))
assert(MKF2.Ts == Ts)
assert(isequaln(MKF2.u_meas, u_meas))
assert(isequal(size(MKF2.Q), [1 4]))
assert(isequal(MKF2.Q{1}, ...
    diag([0.01 0.01 sigma_wp(1, 1)^2/MKF2.d sigma_wp(2, 1)^2/MKF2.d])))
assert(isequal(MKF2.Q{2}, ...
    diag([0.01 0.01 sigma_wp(1, 1)^2/MKF2.d sigma_wp(2, 2)^2/MKF2.d])))
assert(isequal(MKF2.Q{3}, ...
    diag([0.01 0.01 sigma_wp(1, 2)^2/MKF2.d sigma_wp(2, 1)^2/MKF2.d])))
assert(isequal(MKF2.Q{4}, ...
    diag([0.01 0.01 sigma_wp(1, 2)^2/MKF2.d sigma_wp(2, 2)^2/MKF2.d])))
assert(isequal(size(MKF2.R), [1 4]))
assert(isequal(MKF2.R{1}, R) && isequal(MKF2.R{2}, R))
assert(isequal(MKF2.R{3}, R) && isequal(MKF2.R{4}, R))
assert(numel(MKF2.filters) == MKF2.n_filt)
assert(isequal(size(MKF2.seq), [MKF2.n_filt 1]))
assert(isequal(size(cell2mat(MKF2.seq)), [MKF2.n_filt MKF2.f]))
assert(MKF2.beta == sum(MKF2.p_seq))
assert(MKF2.f == size(MKF2.seq{1}, 2))
assert(isequal(size(MKF2.xkp1_est), [n 1]))
assert(isequal(size(MKF2.ykp1_est), [ny 1]))
assert(sum(MKF2.p_gamma) == 1)
alpha = (1 - (1 - epsilon).^d);  % prob. over detection interval 
p_gamma = [1-alpha'; alpha'];
Z = [0 0; 0 1; 1 0; 1 1];  % combinations
p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
p_gamma = p_gamma ./ sum(p_gamma);  % normalized
assert(isequal(round(p_gamma, 6), ...
    [0.960596; 0.019504; 0.019504; 0.000396]))
assert(isequal(round(MKF2.p_gamma, 6), ...
    [0.960596; 0.019504; 0.019504; 0.000396]))

% Check optional definition with an initial state estimate works
x0 = [0.1; 0.5; -0.2; -0.4];
MKF_testx0 = mkf_observer_RODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label,x0);
assert(isequal(MKF_testx0.xkp1_est, x0))
assert(isequal(MKF_testx0.ykp1_est, C * x0))


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

% Custom MKF test observer

% Devise a custom multi-model filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% this test simulation (t = t_shock)
% Multiple model filter 1
A2 = repmat({A}, 1, 3);
Bu2 = repmat({Bu}, 1, 3);
C2 = repmat({C}, 1, 3);
Du2 = repmat({Du}, 1, 3);
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 1 1]);
P0_init = repmat({P0}, 1, 3);
Q2 = {diag([Q0(1,1) Q0(2,2) sigma_wp(1,1)^2 sigma_wp(2,1)^2]), ...
      diag([Q0(1,1) Q0(2,2) sigma_wp(1,2)^2 sigma_wp(2,1)^2]), ...
      diag([Q0(1,1) Q0(2,2) sigma_wp(1,1)^2 sigma_wp(2,2)^2])};
R2 = {diag(sigma_M.^2), diag(sigma_M.^2), diag(sigma_M.^2)};
seq = {zeros(1, nT+1); zeros(1, nT+1); zeros(1, nT+1)};
seq{2}(t == t_shock(1)) = 1;
seq{3}(t == t_shock(2)) = 2;
p_gamma = [1-epsilon epsilon]';
Z = [0 0; 0 1; 1 0];  % combinations
p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
p_gamma = p_gamma ./ sum(p_gamma);  % normalized
T = repmat(p_gamma', 3, 1);
d = 1;
MKF3 = mkf_observer(A2,Bu2,C2,Du2,Ts,P0_init,Q2,R2,seq,T,d,'MKF3');

seq = {zeros(1, nT+1)};
seq{1}(t == t_shock(1)) = 1;
seq{1}(t == t_shock(2)) = 2;
MKF4 = mkf_observer(A2,Bu2,C2,Du2,Ts,P0_init,Q2,R2,seq,T,d,'MKF4');


% Choose observers to test
observers = {KF3, SKF, MKF1, MKF2, MKF3, MKF4};


% Note: KF1 is too slow to pass static error test here

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

% Simulate observers

% Measured inputs (not including disturbances)
U_m = U;

n_obs = numel(observers);
MSE = containers.Map();
for i = 1:n_obs

    obs = observers{i};
    [obs, sim_results] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m,obs, ...
        alpha,u_meas);

    % Check observer errors are zero prior to
    % input disturbance
    assert(all(abs(sim_results.X_est(1:5,:) - X(1:5, :)) < 1e-10, [1 2]))
    assert(all(abs(sim_results.Y_est(1:5,:) - Y(1:5, :)) < 1e-10, [1 2]))

    % Check observer static errors are small
    % after input disturbance
    % TODO: Should these be closer?
    if all(sigma_MP == 0)
        assert(all(abs(sim_results.Y_est(end, :) - Y(end, :)) < 1e-3, [1 2]));
        assert(all(abs(sim_results.X_est(end, 3:4) - du0) < 1e-3, [1 2]));
    end

    % Compute mean-squared error
    Y_est = sim_results.Y_est;
    MSE(obs.label) = mean((Y_est - Y).^2);
    %fprintf("%d, %s: %f\n", i, obs.label, mean((Y_est - Y).^2))

    % Save updated observer
    observers{i} = obs;

end


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
 {'KF3',               'MKF1',              'MKF2',              ...
  'MKF3',              'MKF4',              'SKF'}, ...
 {[0.000676 0.000936], [0.001826 0.006518], [0.002042 0.003289], ...
  [0.001538 0.001718], [0.000123 0.000132], [0.000123 0.000132]} ...
);

for label = MSE.keys
   assert(isequal(round(MSE(label{1}), 6), MSE_test_values(label{1})))
end

return


function [obs, sim_results] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m, ...
    obs,alpha,u_meas)

    k = (0:nT)';
    t = Ts*k;
    X_est = nan(nT+1,n);
    Y_est = nan(nT+1,ny);
    E_obs = nan(nT+1,ny);
    
    % Arrays to store observer variables
    switch obs.label
        case {'MKF1', 'MKF2'}
            n_filt = obs.n_filt;
            MKF_p_seq_g_Yk = nan(nT+1, n_filt);
        otherwise
            n_filt = 1;
    end
    K_obs = cell(nT+1, n_filt);
    trP_obs = cell(nT+1, n_filt);

    % Start simulation at k = 0
    for i = 1:nT+1

        % For debugging:
        %fprintf("t = %f\n", t(i));

        % Process measurements
        uk_m = U_m(i,:)';
        yk_m = Y_m(i,:)';

        % Record observer estimates and output errors
        X_est(i, :) = obs.xkp1_est';
        Y_est(i, :) = obs.ykp1_est';
        E_obs(i, :) = yk_m' - obs.ykp1_est';

        % Kalman update equations
        % Update observer gains and covariance matrix
        switch obs.label

            case {'KF1', 'KF2', 'KF3'}
                obs = update_KF(obs, uk_m, yk_m);

                % Record filter gain and covariance matrix
                K_obs{i, 1} = obs.K(:)';
                trP_obs{i, 1} = trace(obs.P);

            case {'SKF'}

                % Get actual shock occurence indicators
                alpha_k = alpha(i, :);

                % Update observer estimates
                obs = update_SKF(obs, uk_m, yk_m, alpha_k);

                % Record filter gain and covariance matrix
                K_obs{i, 1} = obs.K(:)';
                trP_obs{i, 1} = trace(obs.P);

            case {'MKF1', 'MKF2', 'MKF3', 'MKF4'}
                obs = update_MKF(obs, uk_m, yk_m);

                % Record filter gains and covariance matrices
                for j = 1:obs.n_filt
                    K_obs{i, j} = obs.filters{j}.K(:)';
                    trP_obs{i, j} = trace(obs.filters{j}.P);
                end

                % Record filter conditional probabilities
                MKF_p_seq_g_Yk(i, :) = obs.p_seq_g_Yk';

            otherwise
                error("Value error: observer not recognized")

        end

    end

    sim_results.t = t;
    sim_results.k = k;
    sim_results.X_est = X_est;
    sim_results.Y_est = Y_est;
    sim_results.E_obs = E_obs;
    sim_results.K_obs = K_obs;
    sim_results.trP_obs = trP_obs;
    switch obs.label
        case {'MKF1', 'MKF2'}
            sim_results.MKF_p_seq_g_Yk = MKF_p_seq_g_Yk;
    end

end