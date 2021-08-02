% Test functions mkf_filter_RODD.m and update_MKF.m

clear all
plot_dir = 'plots';

seed = 0;
rng(seed)

% Load system and disturbance model from file
sys_rodin_step

% Set noise variances for observer design
sigma_M = 0.1;
sigma_W = [0; 0];

% Load observers from file
obs_rodin_step

assert(MKF1.epsilon == 0.01)
assert(isequal(MKF1.sigma_wp, sigma_wp))
assert(MKF1.n_filt == 7)
assert(MKF1.i == 0)
assert(MKF1.n == n)
assert(MKF1.nu == nu)
assert(MKF1.ny == ny)
assert(MKF1.nj == 2)
assert(isequal(MKF1.A{1}, A) && isequal(MKF1.A{2}, A))
assert(isequal(MKF1.B{1}, B) && isequal(MKF1.B{2}, B))
assert(isequal(MKF1.C{1}, C) && isequal(MKF1.C{2}, C))
assert(isequal(MKF1.D{1}, D) && isequal(MKF1.D{2}, D))
assert(MKF1.Ts == Ts)
assert(isequal(MKF1.Q{1}, [0.01 0; 0 sigma_wp(1)^2]))
assert(isequal(MKF1.Q{2}, [0.01 0; 0 sigma_wp(2)^2]))
assert(isequal(MKF1.R{1}, R) && isequal(MKF1.R{2}, R))
assert(numel(MKF1.filters) == MKF1.n_filt)
assert(isequal(size(MKF1.seq), [MKF1.n_filt 1]))
assert(isequal(size(cell2mat(MKF1.seq)), [MKF1.n_filt MKF1.f*MKF1.d]))
assert(MKF1.beta == sum(MKF1.p_seq))
assert(MKF1.nf == size(MKF1.seq{1}, 2))
assert(isequal(size(MKF1.xkp1_est), [n 1]))
assert(isequal(size(MKF1.ykp1_est), [ny 1]))
assert(isequal(MKF1.p_gamma, [1-MKF1.epsilon; MKF1.epsilon]))

assert(MKF2.epsilon == 0.01)
assert(MKF2.n_filt == 16)
assert(MKF2.i == 0)
assert(MKF2.n == n)
assert(MKF2.nu == nu)
assert(MKF2.ny == ny)
assert(MKF2.nj == 2)
assert(isequal(MKF2.A{1}, A) && isequal(MKF2.A{2}, A))
assert(isequal(MKF2.B{1}, B) && isequal(MKF2.B{2}, B))
assert(isequal(MKF2.C{1}, C) && isequal(MKF2.C{2}, C))
assert(isequal(MKF2.D{1}, D) && isequal(MKF2.D{2}, D))
assert(MKF2.Ts == Ts)
assert(isequal(MKF2.Q{1}, [0.01 0; 0 sigma_wp(1)^2]))
assert(isequal(MKF2.Q{2}, [0.01 0; 0 sigma_wp(2)^2]))
assert(isequal(MKF2.R{1}, R) && isequal(MKF2.R{2}, R))
assert(numel(MKF2.filters) == MKF2.n_filt)
assert(isequal(size(MKF2.seq), [MKF2.n_filt 1]))
assert(isequal(size(cell2mat(MKF2.seq)), [MKF2.n_filt MKF2.f*MKF2.d]))
assert(MKF2.beta == sum(MKF2.p_seq))
assert(MKF2.nf == size(MKF2.seq{1}, 2))
assert(isequal(size(MKF2.xkp1_est), [n 1]))
assert(isequal(size(MKF2.ykp1_est), [ny 1]))
assert(isequal(MKF2.p_gamma, [1-MKF2.epsilon; MKF2.epsilon]))


% Simulation settings
nT = 100;
t = Ts*(0:nT)';

% Choose time and amplitude of input disturbance
t_shock = 10;
du0 = 1;

% Inputs
%U = (idinput(size(t)) + 1)/2;
U = zeros(size(t));
U(t >= 1) = -1;
alpha = zeros(size(t));
alpha(t == 9.5) = 1;  % this is used by the SKF observer
%Wp = 1*alpha;
Wp = zeros(size(t));  % Set RODD disturbance to 0 for this test
U_sim = [U Wp];

% Apply the input disturbance
Du = zeros(size(U_sim));
Du(t >= t_shock, 1) = du0;


% Custom MKF test observer

% Devise a custom multi-model filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% this test simulation (t = t_shock)
% Multiple model filter 1
A2 = repmat({A}, 1, 2);
B2 = repmat({B}, 1, 2);
C2 = repmat({C}, 1, 2);
D2 = repmat({D}, 1, 2);
P0 = 1000*eye(n);
Q0 = diag([Q1 1]);
P0_init = repmat({P0}, 1, 2);
Q2 = {diag([Q0(1,1) sigma_wp(1,1)^2]), ...
      diag([Q0(1,1) sigma_wp(1,2)^2])};
R2 = {sigma_M^2, sigma_M^2};
seq = {zeros(1, nT+1); zeros(1, nT+1)};
seq{2}(t == 9.5) = 1;
p_gamma = [1-epsilon epsilon]';
MKF3 = mkf_filter(A2,B2,C2,D2,Ts,P0_init,Q2,R2,seq,p_gamma,'MKF3');

seq = {zeros(1, nT+1)};
seq{1}(t == 9.5) = 1;
p_gamma = [1-epsilon epsilon]';
MKF4 = mkf_filter(A2,B2,C2,D2,Ts,P0_init,Q2,R2,seq,p_gamma,'MKF4');

% Choose observers to test
observers = {KF2, KF3, SKF, MKF1, MKF2, MKF3, MKF4};

% Note: KF1 is too slow to pass static error test here


% Simulate system
X = zeros(nT+1,n);
Y = zeros(nT+1,ny);
xk = zeros(n,1);

for i=1:nT+1

    % Inputs
    uk = U_sim(i,:)' + Du(i,:)';

    % Compute y(k)
    yk = C*xk + D*uk;

    % Store results
    X(i, :) = xk';
    Y(i, :) = yk';
    
    % Compute x(k+1)
    xk = A*xk + B*uk;

end

% Check simulation output is correct
[Y2, t, X2] = lsim(Gpss,U_sim + Du,t);
assert(isequal(X, X2))
assert(isequal(Y, Y2))

% Choose measurement noise for plant
sigma_MP = 0;
Y_m = Y + sigma_MP'.*randn(size(Y));


% Simulate observers

% Set disturbance inputs to zero for observer
% simulation since they are not measured.
U_obs = [U zeros(nT+1,1)];

n_obs = numel(observers);
MSE = containers.Map();

for i = 1:n_obs

    obs = observers{i};
    [obs, sim_results] = run_test_simulation(nT,Ts,n,ny,U_obs,Y_m,obs,alpha);

    % Check observer errors are zero prior to
    % input disturbance
    assert(all(abs(sim_results.X_est(1:20,:) - X(1:20, :)) < 1e-10, [1 2]))
    assert(all(abs(sim_results.Y_est(1:20,:) - Y(1:20, :)) < 1e-10))

    % Check observer static errors are small
    % after input disturbance
    if all(sigma_MP == 0)
        assert(abs(sim_results.Y_est(end, :) - Y(end, :)) < 1e-4);
        assert(abs(sim_results.X_est(end, 2) - du0) < 1e-4);
    end
    
    % Compute mean-squared error
    Y_est = sim_results.Y_est;
    MSE(obs.label) = mean((Y_est - Y).^2);
    %fprintf("Observer %s tested\n", obs.label)

end

MSE_test_values = containers.Map(...
    {'KF2', 'KF3', 'SKF', 'MKF1', 'MKF2', 'MKF3', 'MKF4'}, ...
    [0.000934 0.003524 0.000929 0.004914 0.005950 0.002709 0.000929]' ...
);

for label = MSE.keys
   assert(isequal(round(MSE(label{1}), 6), MSE_test_values(label{1})))
end


% % Display results of last simulation
% 
% X_est = sim_results.X_est;
% E_obs = sim_results.E_obs;
% K_obs = sim_results.K_obs;
% trP_obs = sim_results.trP_obs;
% 
% table(t,alpha,U,Du,Wp,X,Y,Y_m,X_est,Y_est,E_obs)
% 
% % Display gains and trace of covariance matrix
% table(t, cell2mat(K_obs), cell2mat(trP_obs), ...
%     'VariableNames',{'t', 'K{1}, K{2}', 'trace(P{1}), trace(P{2})'})
% 
% % Show table of mean-squared errors
% table(MSE.keys', cell2mat(MSE.values'), ...
%     'VariableNames', {'Observer', 'MSE'})
% 
% 
% % Plot of inputs and outputs
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
% title('Process output measurements and estimates')
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
% title('Actual states and estimates')
% grid on
% 
% ax3 = subplot(4,1,3);
% stairs(t,U,'Linewidth',2); hold on;
% stairs(t,Wp,'Linewidth',2)
% stairs(t,Du(:,1),'Linewidth',2)
% max_min = [min(min([U Wp Du(:,1)])) max(max([U Wp Du(:,1)]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$u(k)$, $w_p(k)$ and $d_u(k)$')
% legend('$u(k)$', '$w_p(k)$', '$d_u(k)$')
% title('Actual process inputs')
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
%         filename = sprintf('rod_mkf_filter_test_pyk_wfplot.png');
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
%         filename = sprintf('rod_mkf_filter_test_trP_wfplot.png');
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
%         filename = sprintf('rod_mkf_filter_test_K_wfplot.png');
%         filepath = fullfile(plot_dir, filename);
%         show_waterfall_plot(t, K1_obs, [0 5], ax_labels, [0 82], filepath);
%         title('Filter correction gains (k1)')
% 
% end


function [obs, sim_results] = run_test_simulation(nT,Ts,n,ny,U,Y_m, ...
    obs,alpha)

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
        uk = U(i,:)';
        yk = Y_m(i,:)';

        % Record observer estimates and output errors
        X_est(i, :) = obs.xkp1_est';
        Y_est(i, :) = obs.ykp1_est';
        E_obs(i, :) = yk' - obs.ykp1_est';

        % Kalman update equations
        % Update observer gains and covariance matrix
        switch obs.label

            case {'KF1', 'KF2', 'KF3'}
                obs = update_KF(obs, uk, yk);

                % Record filter gain and covariance matrix
                K_obs{i, 1} = obs.K';
                trP_obs{i, 1} = trace(obs.P);

            case {'SKF'}  % Scheduled Kalman filters

                % Set process noise covariance matrix Q based on
                % actual shock occurence
                a = alpha(i, :);
                n_dist = size(a, 2);
                x_var = diag(obs.Q0);
                u_meas = [1; 0];
                x_var(~u_meas) = obs.sigma_wp(sub2ind(size(obs.sigma_wp), ...
                    1:n_dist, a+1)).^2;
                obs.Q = diag(x_var);

                % Update observer gains and covariance matrix
                obs = update_KF(obs, uk, yk);

                % Record filter gain and covariance matrix
                K_obs{i, 1} = obs.K';
                trP_obs{i, 1} = trace(obs.P);

            case {'MKF1', 'MKF2', 'MKF3', 'MKF4'}
                obs = update_MKF(obs, uk, yk);

                % Record filter gains and covariance matrices
                for j=1:obs.n_filt
                    K_obs{i, j} = obs.filters{j}.K';
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