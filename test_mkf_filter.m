% Test functions mkf_filter.m and update_MKF.m

clear all

seed = 0;
rng(seed)

% Sample period
Ts = 0.5;

% Input disturbance variance
%sigma_w = 0.1;
sigma_w = 0;

% Process noise std. dev.
sigma_W = [0; 0];

% Measurement noise std. dev.
sigma_M = 0;

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

assert(isequal(size(A1), size(A2)))
assert(isequal(size(B1), size(B2)))
assert(isequal(size(C1), size(C2)))
assert(isequal(size(D1), size(D2)))


% Simulate switching system

% System models
A = {A1, A2};
B = {B1, B2};
C = {C1, C2};
D = {D1, D2};

% Simulation settings
nT = 60;
t = Ts*(0:nT)';

% Random inputs
%U = (idinput(size(t)) + 1)/2;
U = zeros(nT+1,1);
U(t>2) = 1;
V = sigma_M * randn(size(t));

% Switching sequence
%Gamma = int8(rand(nT+1, 1) > T(1, 1));
Gamma = int8(zeros(nT+1, 1));
Gamma(t>=10, 1) = 1;

% Simulate system
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

% Plot of inputs and outputs
figure(1); clf

ax1 = subplot(5,1,1:2);
plot(t,Y,'Linewidth',2); hold on
plot(t,Ym,'o');
max_min = [min(min([Y Ym])) max(max([Y Ym]))];
bd = max([0.1 diff(max_min)*0.1]);
ylim(max_min + [-bd bd])
xlabel('t')
ylabel('y(k)')
title('System output and output measurements')
grid on

ax2 = subplot(5,1,3:4);
stairs(t,U,'Linewidth',2);
max_min = [min(min(U)) max(max(U))];
bd = max([0.1 diff(max_min)*0.1]);
ylim(max_min + [-bd bd])
xlabel('t')
ylabel('u(k) and w_p(k)')
legend('u(k)')
title('Input')
grid on

ax3 = subplot(5,1,5);
stairs(t,Gamma,'Linewidth',2)
max_min = [min(min(Gamma)) max(max(Gamma))];
bd = max([0.1 diff(max_min)*0.1]);
ylim(max_min + [-bd bd])
xlabel('t')
ylabel('gamma(k)')
title('Model sequence')
grid on

linkaxes([ax1 ax2 ax3], 'x')


% Kalman filters

% Observer parameters
Q1 = 0.01;
R1 = 0.1^2;
Q2 = 0.01;
R2 = 0.1^2;
Q = {Q1, Q2};
R = {R1, R2};
assert(isequal(size(Q1), size(Q2)))
assert(isequal(size(R1), size(R2)))

P0 = 10000;
x0 = 0.2;
KF1 = kalman_filter(A1,B1,C1,D1,Ts,P0,Q1,R1,'KF1',x0)
KF2 = kalman_filter(A2,B2,C2,D2,Ts,P0,Q2,R2,'KF2',x0)


% Observers with switching systems

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
P0 = repmat({10000}, n_filt, 1);
x0 = 1;
d = 1;
MKF1 = mkf_filter(A,B,C,D,Ts,P0,Q,R,seq,T,d,"MKF1",x0);
assert(isequal(MKF1.xkp1_est, x0))
assert(isequal(MKF1.ykp1_est, C{1} * x0))

% Re-define with no initial state specified (should be set to zero)
MKF1 = mkf_filter(A,B,C,D,Ts,P0,Q,R,seq,T,d,"MKF1");

assert(isequal(MKF1.A, A))
assert(isequal(MKF1.B, B))
assert(isequal(MKF1.C, C))
assert(isequal(MKF1.D, D))
assert(isequal(MKF1.Ts, Ts))
assert(isequal(MKF1.Q, Q))
assert(isequal(MKF1.R, R))
assert(isequal(MKF1.seq, seq))
assert(MKF1.n_filt == n_filt)
assert(MKF1.n_filt == numel(MKF1.filters))
assert(isequaln(MKF1.i, nan(1, 2)))
assert(isequal(MKF1.i_next, int16([1 1])))
assert(MKF1.n == n)
assert(MKF1.nu == nu)
assert(MKF1.ny == ny)
assert(MKF1.f == size(MKF1.seq{1}, 2))
assert(MKF1.nj == 2)
assert(isequal(MKF1.T, T))
assert(isequal(MKF1.xkp1_est, zeros(n, 1)))
assert(MKF1.ykp1_est == 0)


% Simulate observers

observers = {KF1, KF2, MKF1};
n_obs = numel(observers);

X_est = zeros(nT+1, n*n_obs);
Y_est = zeros(nT+1, ny*n_obs);
E_obs = zeros(nT+1, ny*n_obs);

obs_mkf = observers{3};
n_filters = size(obs_mkf.seq, 1);
K_obs = cell(nT+1, n*n_filters);
trP_obs = cell(nT+1, n_filters);
MKF_data = nan(nT+1, 2+n_filters);

for i = 1:nT+1

    yk = Y(i,:)';
    if i > 1
        ukm1 = U(i-1, 1);
    else
        ukm1 = 0;
    end

    % Update observers
    x_est = nan(1, n*n_obs);
    y_est = nan(1, ny*n_obs);
    for f = 1:n_obs
        obs = observers{f};
        if strcmp(obs.label, obs_mkf.label)
            show_plots = false;
            obs = update_MKF2(obs, ukm1, yk, show_plots);
            MKF_data(i, :) = [obs.i obs.p_seq_g_Yk'];
            for j = 1:obs.n_filt
                K_obs{i, j} = obs.filters{j}.K';
                trP_obs{i, j} = trace(obs.filters{j}.P);
            end
        else
            switch obs.label
                case {'KF1', 'KF2'}
                    obs = update_KF(obs, ukm1, yk);
                case {'MKF1', 'MKF2', 'MKF3'}
                    obs = update_MKF2(obs, ukm1, yk);
            end
        end
        x_est(1, (f-1)*n+1:f*n) = obs.xkp1_est';
        y_est(1, (f-1)*ny+1:f*ny) = obs.ykp1_est';
        observers{f} = obs;  % save changes
    end

    % Record observer estimates
    X_est(i, :) = x_est;
    Y_est(i, :) = y_est;
    E_obs(i, :) = yk' - y_est;

end

% Combine results
sim_results = table(t,Gamma,U,X,Y,X_est,Y_est,E_obs,K_obs,trP_obs);

% Display results
sim_results(:,{'t', 'Gamma', 'U', 'X', 'X_est', 'Y', 'Y_est'})

sim_results_MKF = [ ...
    table(t) ... 
    array2table(MKF_data(:, 1:2), 'VariableNames', {'i1', 'i2'}) ...
    array2table(MKF_data(:, 3:end)) ...
]

obs_labels = cell(1, n_obs);
for f  = 1:n_obs
    obs_labels{f} = observers{f}.label;
end

figure(2); clf
axs = nan(1, 1+n);

axs(1) = subplot(1+n,1,1);
plot(t,Y_est,'Linewidth',2); hold on
plot(t,Y,'k--');
max_min = [min(min([Y Y_est])) max(max([Y Y_est]))];
bd = max([0.1 diff(max_min)*0.1]);
ylim(max_min + [-bd bd])
xlabel('t')
ylabel('y(k)')
title('Output measurements')
legend([obs_labels {'true'}],'Interpreter','none')
grid on

for i = 1:n
    axs(1+i) = subplot(1+n,1,1+i);
    data = [X(:,i) X_est(:, i:n:(n*n_obs+i-1))];
    max_min = [min(min(data)) max(max(data))];
    bd = max([0.1 diff(max_min)*0.1]);
    plot(t, X_est(:, i:n:(n*n_obs+i-1)),'Linewidth',2); hold on
    plot(t, X(:,i),'k--'); 
    ylim(max_min + [-bd bd])
    xlabel('t')
    y_label = sprintf('x_%d(k)',i);
    ylabel(y_label)
    legend([obs_labels {'true'}],'Interpreter','none')
    title(strcat('States - ', y_label))
    grid on
end

linkaxes(axs, 'x')

% Compare simulation results to saved data
results_dir = 'results';
filename = 'mkf_test_sim_results.csv';
test_sim_results = readtable(fullfile(results_dir, filename));
%assert(all(abs(sim_results{:, 'X_est'} ...
%    - test_sim_results{:, {'X_est_1', 'X_est_2'}}) < 1e-10, [1 2]))

% Compute mean-squared error
mse = mean((Y_est - Y).^2)
assert(isequal(round(mse, 4), [5.0189    0.3907    0.0532]))

