% Test MKFObserverSF class

clear all

% This is only needed for plotting
addpath("~/ml-plot-utils/")
plot_dir = 'plots';
if ~isfolder(plot_dir)
    mkdir(plot_dir)
end

seed = 0;
rng(seed)


%% Test sequence generation 1

% Load SISO system and disturbance model from file
sys_rodin_step

% Define sequence fusion observer
% This example is used in methods section of thesis report
P0 = eye(n);
Q0 = diag([0.01 0]);
R = sigma_M^2;
f = 9;  % fusion horizon
m = 1;  % maximum number of shocks
d = 3;  % spacing parameter
MKF_SF95 = MKFObserverSF95(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,"MKF_SF95");

seq_test = { ...
    [0 0 0 0 0 0 0 0 0]; ...
    [1 0 0 0 0 0 0 0 0]; ...
    [0 0 0 1 0 0 0 0 0]; ...
    [0 0 0 0 0 0 1 0 0] ...
};
assert(isequal(MKF_SF95.seq, seq_test))


%% Test sequence generation 2

% Load SISO system and disturbance model from file
sys_rodin_step

% Define sequence fusion observer
% This example is used in methods section of thesis report
P0 = eye(n);
Q0 = diag([0.01 0]);
R = sigma_M^2;
f = 10;  % fusion horizon
m = 2;  % maximum number of shocks
d = 5;  % spacing parameter
MKF_SF95 = MKFObserverSF95(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,"MKF_SF95");

seq_test = { ...
    [0 0 0 0 0 0 0 0 0 0]; ...
    [1 0 0 0 0 0 0 0 0 0]; ...
    [0 0 0 0 0 1 0 0 0 0]; ...
    [1 0 0 0 0 1 0 0 0 0] ...
};
assert(isequal(MKF_SF95.seq, seq_test))


%% Test observer initialization - SISO system

% Load SISO system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

test_seq = { ...
    [0 0 0]; ...
    [1 0 0]; ...
    [0 1 0]; ...
    [0 0 1] ...
};

% Check observer attributes
assert(strcmp(MKF_SF1.type, "MKF_SF"))
assert(MKF_SF1.epsilon == epsilon)
assert(isequal(MKF_SF1.sigma_wp, sigma_wp))
assert(MKF_SF1.n_filt == 4)
assert(isequaln(MKF_SF1.i, [0 0]))
assert(MKF_SF1.n == 2)
assert(MKF_SF1.nu == 1)
assert(MKF_SF1.ny == 1)
assert(MKF_SF1.nj == 2)
assert(isequal(MKF_SF1.A{1}, A) && isequal(MKF_SF1.A{2}, A))
assert(isequal(MKF_SF1.B{1}, Bu) && isequal(MKF_SF1.B{2}, Bu))
assert(isequal(MKF_SF1.C{1}, C) && isequal(MKF_SF1.C{2}, C))
assert(isequal(MKF_SF1.D{1}, Du) && isequal(MKF_SF1.D{2}, Du))
assert(MKF_SF1.Ts == Ts)
assert(isequaln(MKF_SF1.u_meas, u_meas))
assert(isequal(size(MKF_SF1.Q), [1 2]))
assert(isequal(MKF_SF1.Q{1}, [0.01 0; 0 sigma_wp(1)^2/MKF_SF1.d]))
assert(isequal(MKF_SF1.Q{2}, [0.01 0; 0 sigma_wp(2)^2/MKF_SF1.d]))
assert(isequal(size(MKF_SF1.R), [1 2]))
assert(isequal(MKF_SF1.R{1}, R) && isequal(MKF_SF1.R{2}, R))
assert(numel(MKF_SF1.filters) == MKF_SF1.n_filt)
assert(isequal(size(MKF_SF1.seq), [MKF_SF1.n_filt 1]))
assert(isequal(size(cell2mat(MKF_SF1.seq)), [MKF_SF1.n_filt MKF_SF1.n_di]))
assert(isequal(MKF_SF1.seq, test_seq))
assert(MKF_SF1.beta == sum(MKF_SF1.p_seq))
assert(MKF_SF1.f == MKF_SF1.d * MKF_SF1.n_di)
assert(isequal(MKF_SF1.xkp1_est, zeros(n,1)))
assert(isequal(MKF_SF1.P, 1000*eye(2)))
assert(isequal(MKF_SF1.ykp1_est, zeros(ny,1)))
alpha = (1 - (1 - MKF_SF1.epsilon).^MKF_SF1.d);  % prob. over detection interval 
p_gamma = [1-alpha'; alpha'];
assert(isequal(round(alpha, 4), 0.0490))
assert(isequal(round(p_gamma, 4), [0.9510; 0.0490]))
assert(isequal(round(MKF_SF1.alpha, 4), round(alpha, 4)))
assert(isequal(round(MKF_SF1.p_gamma, 4), round(p_gamma, 4)))

test_seq = { ...
    [0 0 0 0 0]; ...
    [1 0 0 0 0]; ...
    [0 1 0 0 0]; ...
    [0 0 1 0 0]; ...
    [0 0 0 1 0]; ...
    [0 0 0 0 1]; ...
    [1 1 0 0 0]; ...
    [1 0 1 0 0]; ...
    [1 0 0 1 0]; ...
    [1 0 0 0 1]; ...
    [0 1 1 0 0]; ...
    [0 1 0 1 0]; ...
    [0 1 0 0 1]; ...
    [0 0 1 1 0]; ...
    [0 0 1 0 1]; ...
    [0 0 0 1 1] ...
};

assert(strcmp(MKF_SF2.type, "MKF_SF"))
assert(MKF_SF2.epsilon == epsilon)
assert(isequal(MKF_SF2.sigma_wp, sigma_wp))
assert(MKF_SF2.n_filt == 16)
assert(isequaln(MKF_SF1.i, [0 0]))
assert(MKF_SF2.n == 2)
assert(MKF_SF2.nu == 1)
assert(MKF_SF2.ny == 1)
assert(MKF_SF2.nj == 2)
assert(isequal(MKF_SF2.A{1}, A) && isequal(MKF_SF2.A{2}, A))
assert(isequal(MKF_SF2.B{1}, Bu) && isequal(MKF_SF2.B{2}, Bu))
assert(isequal(MKF_SF2.C{1}, C) && isequal(MKF_SF2.C{2}, C))
assert(isequal(MKF_SF2.D{1}, Du) && isequal(MKF_SF2.D{2}, Du))
assert(MKF_SF2.Ts == Ts)
assert(isequaln(MKF_SF2.u_meas, u_meas))
assert(isequal(size(MKF_SF2.Q), [1 2]))
assert(isequal(MKF_SF2.Q{1}, [0.01 0; 0 sigma_wp(1)^2/MKF_SF2.d]))
assert(isequal(MKF_SF2.Q{2}, [0.01 0; 0 sigma_wp(2)^2/MKF_SF2.d]))
assert(isequal(size(MKF_SF1.R), [1 2]))
assert(isequal(MKF_SF2.R{1}, R) && isequal(MKF_SF2.R{2}, R))
assert(numel(MKF_SF2.filters) == MKF_SF2.n_filt)
assert(isequal(size(MKF_SF2.seq), [MKF_SF2.n_filt 1]))
assert(isequal(size(cell2mat(MKF_SF2.seq)), [MKF_SF2.n_filt MKF_SF2.n_di]))
assert(isequal(MKF_SF2.seq, test_seq))
assert(MKF_SF2.beta == sum(MKF_SF2.p_seq))
assert(MKF_SF2.f == MKF_SF2.d * MKF_SF2.n_di)
assert(isequal(MKF_SF2.xkp1_est, zeros(n,1)))
assert(isequal(MKF_SF2.P, 1000*eye(2)))
assert(isequal(MKF_SF2.ykp1_est, zeros(ny,1)))
alpha = (1 - (1 - MKF_SF2.epsilon).^MKF_SF2.d);  % prob. over detection interval 
p_gamma = [1-alpha'; alpha'];
assert(isequal(round(MKF_SF2.alpha, 4), round(alpha, 4)))
assert(isequal(round(MKF_SF2.p_gamma, 4), round(p_gamma, 4)))

% Check optional definition with an initial state estimate
label = 'MKF_testx0';
x0 = [0.1; 0.5];
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
f = 15;  % fusion horizon
m = 2;  % maximum number of shocks
d = 3;  % spacing parameter
MKF_testx0 = MKFObserverSF(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label,x0);
assert(isequal(MKF_testx0.xkp1_est, x0))
assert(isequal(MKF_testx0.ykp1_est, C * x0))


%% Test observer on SISO system with 1 shock

% Load SISO system and disturbance model from file
sys_rodin_step

% Multiple model observer with sequence fusion
P0 = eye(n);
Q0 = diag([0.01 0]);
R = sigma_M^2;
f = 15;  % fusion horizon
m = 1;  % maximum number of shocks
d = 5;  % spacing parameter
label = 'MKF_SF98';
MKF_SF98 = MKFObserverSF(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% MKF_SF95 with same parameters as MKF_SF98
label = 'MKF_SF95';
MKF_SF95 = MKFObserverSF95(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

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

% Custom MKF test observer
% Devise a custom multi-model filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% this test simulation (t == t_shock)

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Du = D(:, u_meas);

% Multiple model filter - two sequences, one empty, one correct
A2 = repmat({A}, 1, 2);
Bu2 = repmat({Bu}, 1, 2);
C2 = repmat({C}, 1, 2);
Du2 = repmat({Du}, 1, 2);
Q0 = diag([0.01 1]);
%P0_init = repmat({P0}, 1, 2);
Q2 = {diag([Q0(1,1) sigma_wp(1,1)^2]), ...
      diag([Q0(1,1) sigma_wp(1,2)^2])};
R2 = {sigma_M.^2, sigma_M.^2};
seq = {zeros(1, nT+1); zeros(1, nT+1)};
seq{2}(t == t_shock) = 1;
p_gamma = [1-epsilon epsilon]';
T = repmat(p_gamma', 2, 1);
MKF3 = MKFObserver(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq,T,'MKF3');

% Multiple model filter - one sequence with correct shock
seq = {zeros(1, nT+1)};
seq{1}(t == t_shock) = 1;
p_gamma = [1-epsilon epsilon]';
T = repmat(p_gamma', 2, 1);
MKF4 = MKFObserver(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq,T,'MKF4');

% Define scheduled Kalman filter
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
% seq = sum(alpha .* 2.^(1:-1:0), 2)';
SKF = MKFObserverSched(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq{1},"SKF");

% Choose observers to test
observers = {MKF3, MKF4, SKF, MKF_SF95, MKF_SF98};

% Note: KF1 is too slow to pass static error test here

% Combine all input signals
U_sim = [U Wp];

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

% Simulate observers

% Choose measurement noise for plant simulation
sigma_MP = 0;  % Set to zero for testing
Y_m = Y + sigma_MP'.*randn(size(Y));

% Measured inputs (not including disturbances)
U_m = U;

n_obs = numel(observers);
sim_results = struct();
RMSE_results = struct();
for i = 1:n_obs

    obs = observers{i};
    [obs, sim_result] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m,obs);

    % Check observer errors are zero prior to
    % input disturbance
    %assert(all(abs(sim_result.X_est(1:20,:) - X(1:20, :)) < 1e-10, [1 2]))
    %assert(all(abs(sim_result.Y_est(1:20,:) - Y(1:20, :)) < 1e-10))

    % Check observer static errors are small
    % after input disturbance
    % TODO: Should these be closer?
    if all(sigma_MP == 0)
        assert(abs(sim_result.Y_est(end, :) - Y(end, :)) < 1e-3);
        assert(abs(sim_result.X_est(end, 2) - du0) < 1e-3);
    end

    % Compute mean-squared error
    Y_est = sim_result.Y_est;
    RMSE_results.(obs.label) = mean((Y_est - Y).^2);
    %fprintf("%d, %s: %f\n", i, obs.label, mean((Y_est - Y).^2))

    % Save results
    sim_results.(obs.label) = sim_result;

    % Save updated observer (not needed if using handle objects)
    observers{i} = obs;

end

% Select simulation results to display and plot
obs_label = "MKF_SF98";
obs_i = find(cellfun(@(x) x.label, observers) == obs_label);
obs = observers{obs_i};
sim_result = sim_results.(obs_label);

X_est = sim_result.X_est;
E_obs = sim_result.E_obs;
trP_obs = sim_result.trP_obs;
trP_obs_j = sim_result.trP_obs_j;
K_obs_j = sim_result.K_obs_j;
table(t,alpha,U,Wp,X,Y,Y_m,X_est,Y_est,E_obs,trP_obs)

% Display gains and trace of covariance matrices of each filter
table(t, cell2mat(K_obs_j), cell2mat(trP_obs_j), ...
    'VariableNames',{'t', 'K{1}, K{2}', 'trace(P{1}), trace(P{2})'});

% Show table of mean-squared errors
disp(RMSE_results)

% Plot of inputs and outputs

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure(1); clf
colors = get(gca,'colororder');
ax1 = subplot(4,1,1);
stairs(t,Y_m); hold on
stairs(t,Y_est,'Linewidth',2);
ax1.ColorOrder = colors(1:size(Y_m,2),:);
max_min = [min(min([Y_m Y_est])) max(max([Y_m Y_est]))];
bd = max([0.1 diff(max_min)*0.1]);
ylim(max_min + [-bd bd])
ylabel('$y_M(k),\hat{y}(k)$')
title_text = 'Outputs and output estimates';
title(escape_latex_chars(title_text))
legend('$y_M(k)$','$\hat{y}(k)$')
grid on

ax2 = subplot(4,1,2);
stairs(t,X); hold on
stairs(t,X_est,'Linewidth',2);
ax2.ColorOrder = colors(size(Y,2)+1:size(Y,2)+size(X,2),:);
max_min = [min(min([X X_est])) max(max([X X_est]))];
bd = max([0.1 diff(max_min)*0.1]);
ylim(max_min + [-bd bd])
ylabel('$x_i(k),\hat{x}_i(k)$')
labels = repmat({''}, 1, n*2);
for i=1:n
    labels{i} = sprintf("$x_{%d}(k)$", i);
end
for i=1:n
    labels{i+n} = sprintf("$%s{x}_{%d}(k)$", '\hat', i);
end
legend(labels)
title('States and state estimates')
grid on

ax3 = subplot(4,1,3);
stairs(t,U,'Linewidth',2); hold on;
stairs(t,Wp,'Linewidth',2)
max_min = [min(min([U Wp])) max(max([U Wp]))];
bd = max([0.1 diff(max_min)*0.1]);
ylim(max_min + [-bd bd])
ylabel('$u(k),w_p(k)$')
legend('$u(k)$', '$w_p(k)$')
title('Process inputs')
grid on

ax4 = subplot(4,1,4);
stairs(t,alpha,'Color',colors(end,:),'Linewidth',2)
max_min = [min(min(alpha)) max(max(alpha))];
bd = max([0.1 diff(max_min)*0.1]);
ylim(max_min + [-bd bd])
xlabel('Time ($t$)', 'Interpreter', 'latex')
ylabel('$\gamma(k)$')
title('Random shock sequence')
grid on

linkaxes([ax1, ax2 ax3 ax4], 'x')
set(gcf,'Position',[100 200 448 560]);
filename = sprintf('rod_MKF_SF_test_sim_%s_ioplot', obs_label);
save_fig_to_pdf(fullfile(plot_dir, filename));

% Plot of conditional filter probabilities
%view_angle = [0 82];
view_angle = [0 42];

% Plot of trace of filter covariance matrices

switch obs.type
    case {'MKF_SF', 'MKF_SF95'}

        t = Ts*(0:nT)';
        trP_obs_j = cell2mat(sim_result.trP_obs_j);
        p_seq_g_Yk = sim_result.MKF_p_seq_g_Yk;

        figure(11); clf
        ax_labels = {'$t$', 'Filter $f$', '$\mathrm{Tr}(\mathrm{P}_f(k))$'};
        subplot(2, 1, 1)
        make_waterfall_plot(t, trP_obs_j, [0 2], ax_labels, view_angle);
        title('Trace of error covariance matrices')

        subplot(2, 1, 2)
        ax_labels = {'Time ($t$)', 'Filter $f$', '$\mathrm{Pr}(\Gamma_f(k) \mid \mathbf{Y}(k))$'};
        make_waterfall_plot(t, p_seq_g_Yk, [0 1], ...
            ax_labels, view_angle);
        title('Conditional probabilities of hypotheses')

end

set(gcf,'Position',[660 450 448 360]);
filename = sprintf('rod_MKF_test_sim_%s_prob.png', obs_label);
exportgraphics(gcf,fullfile(plot_dir, filename),'Resolution',600)
%Waterfall plot does not seem to support pdfs
%save_fig_to_pdf(fullfile(plot_dir, filename));

% Plot of filter correction gains (k1)

switch obs.type
    case {'MKF_SF', 'MKF_SF95'}
        K_obs_j = cell2mat(sim_result.K_obs_j);
        % Two gain values
        K1_obs = K_obs_j(:,1:2:end);
        K2_obs = K_obs_j(:,2:2:end);

        figure(12); clf

        subplot(2, 1, 1)
        t = Ts*(0:nT)';
        ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$K_{1,1}(k)$'};
        make_waterfall_plot(t, K1_obs, [0 5], ax_labels, view_angle);
        title('Filter correction gain $K_{1,1}$', 'Interpreter', 'latex')

        subplot(2, 1, 2)
        ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$K_{2,1}(k)$'};
        make_waterfall_plot(t, K2_obs, [0 3], ax_labels, view_angle);
        title('Filter correction gain $K_{2,1}$', 'Interpreter', 'latex')

        set(gcf,'Position',[660 950 448 360]);
        filename = sprintf('rod_MKF_test_sim_%s_K.png', obs_label);
        exportgraphics(gcf,fullfile(plot_dir, filename),'Resolution',600)

end

% Earlier test results (with a shock of amplitude 1)
% MSE_test_values = containers.Map(...
%     {'KF2',   'KF3',   'MKF_SF1',  'MKF_SF2',  'MKF3',  'MKF4',  "SKF"}, ...
%     [0.000934 0.003524 0.004914 0.005016 0.002709 0.000929 0.000929] ...
% );

% Results on Nov 8 before reverting back the Bayesian updating
%MSE_test_values = containers.Map(...
%  {'KF2',   'KF3',   'MKF_SF1',  'MKF_SF2',  'MKF3',  'MKF4',  "SKF"}, ...
%  [0.000934 0.003524 0.009456 0.005016 0.002709 0.000929 0.000929] ...
%);
% Changes since previous results: 
%  - f, m, d parameters for MK1 changed.
%  - Shock probability and variance modified to reflect detection
%    intervals.
%  - Bayesian prob updates only at end of detection intervals
% Note: performance of MKF3 increases if shock amplitude is increased.

% Test results
MSE_test_values = struct( ...
    'KF2', 0.000934, ...
    'KF3', 0.003524, ...
    'MKF_SF98', 0.010453, ...
    'MKF_SF2', 0.009017, ...
    'MKF_SF95', 0.006052, ...
    'MKF3', 0.002709, ...
    'MKF4', 0.000929, ...
    'SKF', 0.000929 ...
);

for label = fieldnames(RMSE_results)'
    fprintf("%s: %f (%f)\n", label{1}, RMSE_results.(label{1}), MSE_test_values.(label{1}))
    assert(isequal(round(RMSE_results.(label{1}), 6), MSE_test_values.(label{1})))
end


%% Test MKF observers on 2x2 system

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
KF1 = KalmanFilter(A,Bu,C,Du,Ts,P0,Q,R,'KF1');

% Kalman filter 2 - tuned to sigma_wp(2)
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([0.01 0.01 sigma_wp(1,2)^2 sigma_wp(2,2)^2]);
R = diag(sigma_M.^2);
KF2 = KalmanFilter(A,Bu,C,Du,Ts,P0,Q,R,'KF2');

% Kalman filter 3 - manually tuned
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([0.01 0.01 0.1^2 0.1^2]);
R = diag(sigma_M.^2);
KF3 = KalmanFilter(A,Bu,C,Du,Ts,P0,Q,R,'KF3');

% Multiple model filter 1
label = 'MKF_SF1';
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
f = 6;  % fusion horizon
m = 1;  % maximum number of shocks
d = 2;  % spacing parameter
MKF_SF1 = MKFObserverSF(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Multiple model filter 2
label = 'MKF_SF2';
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
f = 10;  % fusion horizon
m = 2;  % maximum number of shocks
d = 2;  % spacing parameter
MKF_SF2 = MKFObserverSF(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Check observer initialization
assert(isequal(MKF_SF1.epsilon, epsilon))
assert(isequal(MKF_SF1.sigma_wp, sigma_wp))
assert(MKF_SF1.n_filt == 7)
assert(isequaln(MKF_SF1.i, [0 0]))
assert(MKF_SF1.n == 4)
assert(MKF_SF1.nu == 2)
assert(MKF_SF1.ny == 2)
assert(MKF_SF1.nj == 3)
assert(isequal(MKF_SF1.A{1}, A) && isequal(MKF_SF1.A{2}, A))
assert(isequal(MKF_SF1.B{1}, Bu) && isequal(MKF_SF1.B{2}, Bu))
assert(isequal(MKF_SF1.C{1}, C) && isequal(MKF_SF1.C{2}, C))
assert(isequal(MKF_SF1.D{1}, Du) && isequal(MKF_SF1.D{2}, Du))
assert(MKF_SF1.Ts == Ts)
assert(isequaln(MKF_SF1.u_meas, u_meas))
assert(isequal(size(MKF_SF1.Q), [1 3]))
assert(isequal(MKF_SF1.Q{1}, ...
    diag([0.01 0.01 sigma_wp(1, 1)^2/MKF_SF1.d sigma_wp(2, 1)^2/MKF_SF1.d])))
assert(isequal(MKF_SF1.Q{2}, ...
    diag([0.01 0.01 sigma_wp(1, 1)^2/MKF_SF1.d sigma_wp(2, 2)^2/MKF_SF1.d])))
assert(isequal(MKF_SF1.Q{3}, ...
    diag([0.01 0.01 sigma_wp(1, 2)^2/MKF_SF1.d sigma_wp(2, 1)^2/MKF_SF1.d])))
assert(isequal(size(MKF_SF1.R), [1 3]))
assert(isequal(MKF_SF1.R{1}, R) && isequal(MKF_SF1.R{2}, R) && isequal(MKF_SF1.R{3}, R))
assert(numel(MKF_SF1.filters) == MKF_SF1.n_filt)
assert(isequal(size(MKF_SF1.seq), [MKF_SF1.n_filt 1]))
assert(isequal(size(cell2mat(MKF_SF1.seq)), [MKF_SF1.n_filt MKF_SF1.n_di]))
assert(MKF_SF1.beta == sum(MKF_SF1.p_seq))
assert(MKF_SF1.f == MKF_SF1.n_di * MKF_SF1.d)
assert(isequal(MKF_SF1.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_SF1.P, 1000*eye(4)))
assert(isequal(MKF_SF1.ykp1_est, zeros(ny, 1)))
assert(sum(MKF_SF1.p_gamma) == 1)
alpha = (1 - (1 - epsilon).^d);  % prob. over detection interval 
p_gamma = [1-alpha'; alpha'];
Z = [0 0; 0 1; 1 0];  % combinations
p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
p_gamma = p_gamma ./ sum(p_gamma);  % normalized
assert(isequal(round(p_gamma, 6), [0.960977; 0.019512; 0.019512]))
assert(isequal(round(MKF_SF1.p_gamma, 6), [0.960977; 0.019512; 0.019512]))

% Check observer initialization
assert(isequal(MKF_SF2.epsilon, epsilon))
assert(isequal(MKF_SF2.sigma_wp, sigma_wp))
assert(MKF_SF2.n_filt == 56)
assert(isequaln(MKF_SF2.i, [0 0]))
assert(MKF_SF2.n == 4)
assert(MKF_SF2.nu == 2)
assert(MKF_SF2.ny == 2)
assert(MKF_SF2.nj == 4)
assert(isequal(MKF_SF2.A{1}, A) && isequal(MKF_SF2.A{2}, A))
assert(isequal(MKF_SF2.B{1}, Bu) && isequal(MKF_SF2.B{2}, Bu))
assert(isequal(MKF_SF2.C{1}, C) && isequal(MKF_SF2.C{2}, C))
assert(isequal(MKF_SF2.D{1}, Du) && isequal(MKF_SF2.D{2}, Du))
assert(MKF_SF2.Ts == Ts)
assert(isequaln(MKF_SF2.u_meas, u_meas))
assert(isequal(size(MKF_SF2.Q), [1 4]))
assert(isequal(MKF_SF2.Q{1}, ...
    diag([0.01 0.01 sigma_wp(1, 1)^2/MKF_SF2.d sigma_wp(2, 1)^2/MKF_SF2.d])))
assert(isequal(MKF_SF2.Q{2}, ...
    diag([0.01 0.01 sigma_wp(1, 1)^2/MKF_SF2.d sigma_wp(2, 2)^2/MKF_SF2.d])))
assert(isequal(MKF_SF2.Q{3}, ...
    diag([0.01 0.01 sigma_wp(1, 2)^2/MKF_SF2.d sigma_wp(2, 1)^2/MKF_SF2.d])))
assert(isequal(MKF_SF2.Q{4}, ...
    diag([0.01 0.01 sigma_wp(1, 2)^2/MKF_SF2.d sigma_wp(2, 2)^2/MKF_SF2.d])))
assert(isequal(size(MKF_SF2.R), [1 4]))
assert(isequal(MKF_SF2.R{1}, R) && isequal(MKF_SF2.R{2}, R))
assert(isequal(MKF_SF2.R{3}, R) && isequal(MKF_SF2.R{4}, R))
assert(numel(MKF_SF2.filters) == MKF_SF2.n_filt)
assert(isequal(size(MKF_SF2.seq), [MKF_SF2.n_filt 1]))
assert(isequal(size(cell2mat(MKF_SF2.seq)), [MKF_SF2.n_filt MKF_SF2.n_di]))
assert(MKF_SF2.beta == sum(MKF_SF2.p_seq))
assert(MKF_SF2.f == MKF_SF2.n_di * MKF_SF2.d)
assert(isequal(MKF_SF1.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_SF1.P, 1000*eye(4)))
assert(isequal(MKF_SF1.ykp1_est, zeros(ny, 1)))
assert(sum(MKF_SF2.p_gamma) == 1)
alpha = (1 - (1 - epsilon).^d);  % prob. over detection interval 
p_gamma = [1-alpha'; alpha'];
Z = [0 0; 0 1; 1 0; 1 1];  % combinations
p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
p_gamma = p_gamma ./ sum(p_gamma);  % normalized
assert(isequal(round(p_gamma, 6), ...
    [0.960596; 0.019504; 0.019504; 0.000396]))
assert(isequal(round(MKF_SF2.p_gamma, 6), ...
    [0.960596; 0.019504; 0.019504; 0.000396]))

% Check optional definition with an initial state estimate works
x0 = [0.1; 0.5; -0.2; -0.4];
MKF_testx0 = MKFObserverSF(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
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
%P0_init = repmat({P0}, 1, 3);
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
MKF3 = MKFObserver(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq,T,'MKF3');

seq = {zeros(1, nT+1)};
seq{1}(t == t_shock(1)) = 1;
seq{1}(t == t_shock(2)) = 2;
MKF4 = MKFObserver(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq,T,'MKF4');

% Define scheduled Kalman filter
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
% seq = sum(alpha .* 2.^(1:-1:0), 2)';
SKF = MKFObserverSched(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq{1},"SKF");

% Choose observers to test
observers = {KF3, SKF, MKF_SF1, MKF_SF2, MKF3, MKF4};

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
RMSE_results = containers.Map();
for i = 1:n_obs

    obs = observers{i};
    [obs, sim_result] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m,obs);

    % Check observer errors are zero prior to
    % input disturbance
    assert(all(abs(sim_result.X_est(1:5,:) - X(1:5, :)) < 1e-10, [1 2]))
    assert(all(abs(sim_result.Y_est(1:5,:) - Y(1:5, :)) < 1e-10, [1 2]))

    % Check observer static errors are small
    % after input disturbance
    % TODO: Should these be closer?
    if all(sigma_MP == 0)
        assert(all(abs(sim_result.Y_est(end, :) - Y(end, :)) < 1e-3, [1 2]));
        assert(all(abs(sim_result.X_est(end, 3:4) - du0) < 1e-3, [1 2]));
    end

    % Compute mean-squared error
    Y_est = sim_result.Y_est;
    RMSE_results(obs.label) = mean((Y_est - Y).^2);
    %fprintf("%d, %s: %f\n", i, obs.label, mean((Y_est - Y).^2))

    % Save updated observer
    observers{i} = obs;

end


% Display results of last simulation
% X_est = sim_results.X_est;
% E_obs = sim_results.E_obs;
% trP_obs = sim_results.trP_obs;
% K_obs_j = sim_results.K_obs_j;
% trP_obs_j = sim_results.trP_obs_j;
% table(t,alpha,U,Wp,X,Y,Y_m,X_est,Y_est,E_obs)

% Display gains and trace of covariance matrix
%K_data = cell2mat(cellfun(@(X) X(:)', K_obs, 'UniformOutput', false));
%table(t, K_data, cell2mat(trP_obs), ...
%    'VariableNames', {'t', 'K{1}, K{2}', 'trace(P{1}), trace(P{2})'})

% Show table of mean-squared errors
% table(MSE.keys', cell2mat(MSE.values'), ...
%     'VariableNames', {'Observer', 'MSE'});

% Results on Nov 8 after reverting back the Bayesian updating
MSE_test_values = containers.Map(...
 {'KF3',               'MKF_SF1',           'MKF_SF2', ...
  'MKF3',              'MKF4',              'SKF'}, ...
 {[0.000676 0.000936], [0.001826 0.006518], [0.002042 0.003289], ...
  [0.001538 0.001718], [0.000123 0.000132], [0.000123 0.000132]} ...
);

for label = RMSE_results.keys
    %fprintf("%s: %f, %f (%f, %f)\n", label{1}, MSE(label{1}), MSE_test_values(label{1}))
    assert(isequal(round(RMSE_results(label{1}), 6), MSE_test_values(label{1})))
end
