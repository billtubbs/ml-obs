% Test disturbance model identification from simulated
% system data

clear all;

addpath("~/ml-plot-utils")

seed = 0;
rng(seed);

plot_dir = 'plots';

% Load system and disturbance model from file
sys_rodin_step

% Change measurement noise
sigma_M = 0.5;

% Load observers from file
obs_rodin_step

% Sequenece length
nT = 1000;

% Simulate system
X0 = zeros(n,1);
t = Ts*(0:nT)';
[Wp, Gamma] = sample_random_shocks(nT+1, epsilon, sigma_wp(2), sigma_wp(1));
U = zeros(nT+1,1);
U(t>=5) = 1;
[Y, t, X] = lsim(Gpss,[U Wp],t,X0);
V = sigma_M*randn(nT+1, 1);
Ym = Y + V;  % measurement

observers = {KF1, KF2, KF3, MKF_SP1, MKF_SP2};
obs_labels = cellfun(@(x) x.label, observers, "UniformOutput", true);
f_mkf = find(strcmp("MKF_SP1", obs_labels));

[Xk_est,Yk_est,DiagPk,MKF_vars] = ...
      run_simulation_obs(Ym,U,Gamma,[],observers,f_mkf);

% Calculate estimation errors
E_obs = Y - Yk_est;
rmse = sqrt(sum(E_obs.^2));
RMSE_table = array2table(rmse, 'VariableNames', obs_labels);
disp(RMSE_table)

figure(1); clf
plot_obs_estimates(t,X,Xk_est,Y,Yk_est,obs_labels)

% Check results are as expected
assert(all(RMSE_table{:, 'KF3'} ...
    < RMSE_table{:, {'KF1', 'KF2'}}))
assert(all(RMSE_table{:, "MKF_SP1"} ...
    < RMSE_table{:, {'KF1', 'KF2', 'KF3'}}))

% Run experiments to see effect of different disturbance
% model parameters on estimation errors

% Paremeter values to test
epsilon_values = logspace(-3, -1, 9);
sigma_wp_values = logspace(-1, 1, 9);

% Search grid
n_epsilon = length(epsilon_values);
n_sigma_wp = length(sigma_wp_values);
[idx1, idx2] = ndgrid(1:n_epsilon, 1:n_sigma_wp);
combs = [reshape(idx1, [], 1) reshape(idx2, [], 1)];
n_combs = size(combs, 1);

disp("Running simulations...")
rmse = nan(n_epsilon, n_sigma_wp);
for i = 1:n_combs

    eps_test = epsilon_values(combs(i, 1));
    sig2_test = sigma_wp_values(combs(i, 2));

    % Multiple model observer with sequence pruning #1
    label = 'MKF_SP';
    P0 = 1000*eye(n);
    Q0 = diag([q1 0]);
    R = sigma_M^2;
    nh = 10;  % number of filters
    n_min = 7;  % minimum life of cloned filters
    MKF_SP = MKFObserverSP(model,u_meas,P0,eps_test, ...
        [sigma_wp(1) sig2_test],Q0,R,nh,n_min,label);

    % Run simulation
    observers = {MKF_SP};
    f_mkf = 1;

    [Xk_est,Yk_est,DiagPk,MKF_vars] = ...
        run_simulation_obs(Ym,U,Gamma,[],observers,f_mkf);

    % Calculate estimation errors
    E_obs = Y - Yk_est;
    rmse(combs(i, 1), combs(i, 2)) = sqrt(sum(E_obs.^2));

    fprintf("Sim %d: (%0.2f, %0.2f) %.3f \n", i, epsilon, sigma_wp(2), rmse(i))

end

% Find simulation with lowest errors
i_min = find(rmse == min(rmse, [], [1 2]));

epsilon_best = epsilon_values(combs(i_min, 1));
sigma_wp_best = sigma_wp_values(combs(i_min, 2));
fprintf("Best result: %d\n  Epsilon: %0.2f\n  Sigma_wp: %0.2f\n", ...
    i_min, epsilon_best, sigma_wp_best)


figure(2); clf
x_pts = log(epsilon_values(combs(:, 1)));
y_pts = log(sigma_wp_values(combs(:, 2)));
plot(x_pts, y_pts, '.'); hold on
set(gca, 'TickLabelInterpreter', 'latex')
X = log(epsilon_values(idx1));
Y = log(sigma_wp_values(idx2));
Z = rmse;
[C,h] = contour(X,Y,Z);
clabel(C,h)
plot(log(epsilon_best), log(sigma_wp_best), 'x')
plot(log(epsilon), log(sigma_wp(2)), 'o')
text(log(epsilon) + 0.1, log(sigma_wp(2)) - 0.1, sprintf("$%.2f$", rmse(i_min)), 'Interpreter', 'latex')
xlabel("$\log{\epsilon}$", 'Interpreter', 'latex')
ylabel("$\log{\sigma_{w_p,2}}$", 'Interpreter', 'latex')
grid on
legend("simulation results", "contours of RMSE", "lowest RMSE", "true values", 'Interpreter', 'latex')

filename = "sysid_contplot";
save_fig_to_pdf(fullfile(plot_dir, filename))