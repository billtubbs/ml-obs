% Test disturbance model identification from simulated
% system data

clear all;

addpath("~/ml-plot-utils")

seed = 0;
rng(seed);

plot_dir = 'plots';
results_dir = 'results';

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
Y_m = Y + V;  % measurement

figure(1); clf
make_iodplot(Y,Y_m,t,[U Wp],{'$u(k)$', '$w_p(k)$'})
filename = "sysid_ioplot";
save_fig_to_pdf(fullfile(plot_dir, filename))

observers = {KF1, KF2, KF3, MKF_SP1, MKF_SP2};
obs_labels = cellfun(@(x) x.label, observers, "UniformOutput", true);
f_mkf = find(strcmp("MKF_SP1", obs_labels));

[Xk_est,Yk_est,DiagPk,MKF_vars] = ...
      run_simulation_obs(Y_m,U,Gamma,[],observers,f_mkf);

% Calculate estimation errors
E_obs = Y - Yk_est;
rmses = sqrt(sum(E_obs.^2));
RMSE_table = array2table(rmses, 'VariableNames', obs_labels);
disp(RMSE_table)

figure(2); clf
select = find(strcmp("MKF_SP2", obs_labels));
plot_obs_estimates(t, X, Xk_est(:, (select-1)*n+1:select*n), Y, Yk_est(:, select), escape_latex_chars(obs_labels(select)))
filename = "sysid_obs_est";             
save_fig_to_pdf(fullfile(plot_dir, filename))

% Check results are as expected
assert(all(RMSE_table{:, 'KF3'} ...
    < RMSE_table{:, {'KF1', 'KF2'}}))
assert(all(RMSE_table{:, "MKF_SP1"} ...
    < RMSE_table{:, {'KF1', 'KF2', 'KF3'}}))

% Run experiments to see effect of different disturbance
% model parameters on estimation errors

% Choose simulation case (defined below)
sim_case = 2;

switch sim_case
    case 1
        % Range of paremeter values to test
        epsilon_values = logspace(-3, -1, 9);
        sigma_wp_values = logspace(-1, 1, 9);
    case 2
        % Range of paremeter values to test
        epsilon_values = logspace(-2.25, -1.75, 9);
        sigma_wp_values = logspace(-0.25, 0.25, 9);
end

n_epsilon = length(epsilon_values);
n_sigma_wp = length(sigma_wp_values);

% Search grid
[idx1, idx2] = ndgrid(1:n_epsilon, 1:n_sigma_wp);
combs = [reshape(idx1, [], 1) reshape(idx2, [], 1)];
n_combs = size(combs, 1);

% Results filename
filename = sprintf("rmse_results_SISO_%d.csv", sim_case);

if ~exist(fullfile(results_dir, filename), 'file')

    fprintf("Running %d simulations...\n", n_combs)
    results = array2table(nan(n_combs, 5), 'VariableNames', {'x', 'y', ...
        'epsilon', 'sigma_wp_2', 'RMSE'});
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
            run_simulation_obs(Y_m,U,Gamma,[],observers,f_mkf);

        % Calculate estimation errors
        E_obs = Y - Yk_est;
        rmse = sqrt(sum(E_obs.^2));
        results{i, {'x', 'y', 'epsilon', 'sigma_wp_2', 'RMSE'}} = ...
            [combs(i, :) eps_test sig2_test rmse];
        fprintf("Sim %d: (%0.2f, %0.2f) %.3f \n", i, eps_test, sig2_test, rmse)
    
    end

    % Save results
    writetable(results, fullfile(results_dir, filename));
    fprintf("Results saved to file '%s'\n", filename)

else

    fprintf("Loading simulation results from file '%s'\n", filename)
    results = readtable(fullfile(results_dir, filename));
    assert(size(results, 1) == n_combs);

end

% Make a grid of RMSE values
rmses = nan(n_epsilon, n_sigma_wp);
for i = 1:n_combs
    comb = results{i, {'x', 'y'}};
    eps_test = results{i, 'epsilon'};
    sig2_test = results{i, 'sigma_wp_2'};
    rmses(comb(1), comb(2)) = results{i, 'RMSE'};
end

% Find simulation with lowest errors
i_min = find(rmses == min(rmses, [], [1 2]));

epsilon_best = epsilon_values(combs(i_min, 1));
sigma_wp_best = sigma_wp_values(combs(i_min, 2));
fprintf("Best result: %d\n  Epsilon: %0.3f\n  Sigma_wp: %0.3f\n", ...
    i_min, epsilon_best, sigma_wp_best)

% Make contour plot
figure(2); clf
x_pts = log(epsilon_values(combs(:, 1)));
y_pts = log(sigma_wp_values(combs(:, 2)));
plot(x_pts, y_pts, '.'); hold on
set(gca, 'TickLabelInterpreter', 'latex')
X = log(epsilon_values(idx1));
Y = log(sigma_wp_values(idx2));
Z = rmses;
[C,h] = contour(X,Y,Z);
clabel(C,h)
plot(log(epsilon_best), log(sigma_wp_best), 'x')
plot(log(epsilon), log(sigma_wp(2)), 'o')
text(log(epsilon) + 0.1, log(sigma_wp(2)) - 0.1, sprintf("$%.2f$", rmses(i_min)), 'Interpreter', 'latex')
xlabel("$\log{\epsilon}$", 'Interpreter', 'latex')
ylabel("$\log{\sigma_{w_p,2}}$", 'Interpreter', 'latex')
grid on
legend("simulation results", "contours of RMSE", "lowest RMSE", ...
    "true values", 'Interpreter', 'latex')

filename = sprintf("sysid_contplot_%d", sim_case);
save_fig_to_pdf(fullfile(plot_dir, filename))