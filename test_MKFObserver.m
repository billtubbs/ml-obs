% Test functions mkf_observer.m and update_MKF.m

clear all

seed = 0;
rng(seed)

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


%% Simulation tests

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


% Kalman filters

% Observer parameters (same for all observers)
P0 = 10000;
x0 = 0.5;
Q1 = 0.01;
R1 = 0.1^2;
Q2 = 0.01;
R2 = 0.1^2;
assert(isequal(size(Q1), size(Q2)))
assert(isequal(size(R1), size(R2)))

KF1 = KalmanFilter(A1,B1,C1,D1,Ts,P0,Q1,R1,'KF1',x0);
KF2 = KalmanFilter(A2,B2,C2,D2,Ts,P0,Q2,R2,'KF2',x0);

% Define observer with a switching system
Q = {Q1,Q2};
R = {R1,R2};

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
P0j = repmat({P0}, n_filt, 1);
d = 1;

% First, define with no initial state specified (should be set to zero)
MKF1 = MKFObserver(A,B,C,D,Ts,P0j,Q,R,seq,T,d,'MKF1');

assert(strcmp(MKF1.type, "MKF"))
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

% Redefine this time with initial conditions
MKF1 = MKFObserver(A,B,C,D,Ts,P0j,Q,R,seq,T,d,'MKF1',x0);
assert(isequal(MKF1.xkp1_est, x0))
assert(isequal(MKF1.ykp1_est, C{1} * x0))
gamma0 = 0;
MKF1 = MKFObserver(A,B,C,D,Ts,P0j,Q,R,seq,T,d,'MKF1',x0,gamma0);
assert(isequal(MKF1.xkp1_est, x0))
assert(isequal(MKF1.ykp1_est, C{1} * x0))
assert(isequal(MKF1.gamma_k, zeros(n_filt, 1)))
gamma0 = zeros(n_filt, 1);
gamma0(end) = 1;
MKF1 = MKFObserver(A,B,C,D,Ts,P0j,Q,R,seq,T,d,'MKF1',x0,gamma0);
assert(isequal(MKF1.xkp1_est, x0))
assert(isequal(MKF1.ykp1_est, C{1} * x0))
assert(isequal(MKF1.gamma_k, gamma0))

% With default initial conditions
MKF1 = MKFObserver(A,B,C,D,Ts,P0j,Q,R,seq,T,d,'MKF1');

% Choose observers to include in simulation
observers = {KF1, KF2, MKF1};
n_obs = numel(observers);
obs_labels = cell(1, n_obs);
for f  = 1:n_obs
    obs_labels{f} = observers{f}.label;
end

% Identify which observer to log MKF data for
f_mkf = 3;

% Simulate observers - without measurement noise (Y)
[Xkp1_est,Ykp1_est,MKF_K_obs,MKF_trP_obs,MKF_i,MKF_p_seq_g_Yk,observers1] = ...
    run_simulation_obs(Y,U,observers,f_mkf);

% Move estimates to correct time instants
X_est = [nan(1,n*n_obs); Xkp1_est(1:end-1,:)];
Y_est = [nan(1,ny*n_obs); Ykp1_est(1:end-1,:)];

% Output estimation errors
E_obs = Y - Y_est;

% Combine and display results
%sim_results1 = table(t,Gamma,U,X,Y,Ym,X_est,Y_est,E_obs)

%figure(2); clf
%plot_obs_estimates(t,X,X_est,Y,Y_est,obs_labels)

% Check KF1 was accurate befomsere system switched
assert(nanmean(E_obs(t < 10, 1).^2) < 0.0001)

% Check KF2 was accurate after system switched
assert(nanmean(E_obs(t > 12, 2).^2) < 0.001)

% Compute mean-squared error
mses = nanmean(E_obs.^2);

% Check MKF observer estimation error was low
assert(round(mses(f_mkf), 4) == 0.1296)

% Reset observer states to original initial conditions
MKF1.reset()

% Redefine a new observer (identical to above)
MKF1_new = MKFObserver(A,B,C,D,Ts,P0j,Q,R,seq,T,d,'MKF1');
assert(isequaln(MKF1_new, MKF1))
MKF1_new.label = "MKF1_new";

% Make a copy
MKF1_copy = MKF1_new.copy();
assert(isequaln(MKF1_copy, MKF1_new))
MKF1_copy.label = "MKF1_copy";

% Choose observers to include in simulation
observers = {KF1, KF2, MKF1, MKF1_new, MKF1_copy};

n_obs = numel(observers);
obs_labels = cell(1, n_obs);
for f  = 1:n_obs
    obs_labels{f} = observers{f}.label;
end

% Identify which observer to log MKF data for
f_mkf = 3;

% Simulate observers - with measurement noise (Ym)
[Xkp1_est,Ykp1_est,MKF_K_obs,MKF_trP_obs,MKF_i,MKF_p_seq_g_Yk,observers] = ...
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

% Compute mean-squared error
mses = nanmean(E_obs.^2);
%array2table(mses,'VariableNames',obs_labels)

% Check MKF observer estimation error
assert(round(mses(f_mkf), 4) == 0.1335)

% mses =
%    5.1748    0.6125    0.1335    0.1387    0.0877
% TODO: Why do the copies not produce identical simulation results?
% (see plot figure).


%% Test copy methods

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

% Define MKF observer 1
seq = seq1;
n_filt = numel(seq);
P0j = repmat({P0}, n_filt, 1);
d = 1;

% Define multi-model observer with initial conditions
MKF = MKFObserver(A,B,C,D,Ts,P0j,Q,R,seq,T,d,'MKF1',x0);

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


function [Xkp1_est,Ykp1_est,MKF_K_obs,MKF_trP_obs,MKF_i,MKF_p_seq_g_Yk,observers] = ...
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
    MKF_K_obs = cell(nT+1, n*n_filters);
    MKF_trP_obs = nan(nT+1, n_filters);
    MKF_i = nan(nT+1, 2);
    MKF_p_seq_g_Yk = nan(nT+1, n_filters);

    for i = 1:nT+1

        yk = Ym(i, :)';
        uk = U(i, :);

        % Update observers
        for f = 1:n_obs
            obs = observers{f};
            switch obs.type
                case 'KF'
                    obs.update(yk, uk);
                case 'MKF'
                    obs.update(yk, uk, false);
                otherwise
                    error('Observer type not valid')
            end
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
        end

        % Record observer estimates
        Xkp1_est(i, :) = xkp1_est;
        Ykp1_est(i, :) = ykp1_est;

    end
end


function plot_obs_estimates(t,X,X_est,Y,Y_est,obs_labels)
% Display plots of observer estimates compared to
% true values.
    n_obs = numel(obs_labels);
    n = size(X_est, 2) / n_obs;
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

end