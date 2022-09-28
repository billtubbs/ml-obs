% Test disturbance model identification from simulated
% system data

clear all;

addpath("~/ml-plot-utils")

seed = 0;
rng(seed);

plot_dir = 'plots';
results_dir = 'results';

% Load system and disturbance model from file
sys_rodin_stepramp

% Change measurement noise
sigma_M = 0.5;

% Load observers from file
obs_rodin_stepramp

% Sequenece length
nT = 1000;

% Simulate system
X0 = zeros(n,1);
t = Ts*(0:nT)';
[Wp, Gamma] = sample_random_shocks([nT+1 nw], epsilon, sigma_wp(:, 2), ...
    sigma_wp(:, 1));
U = zeros(nT+1,1);
U(t>=5) = 1;
[Y, t, X] = lsim(Gpss,[U Wp],t,X0);
V = sigma_M*randn(nT+1, 1);
Y_m = Y + V;  % measurement

figure(1); clf
if nu > 1
    u_labels = compose("$u_%d(k)$", 1:nu);
else
    u_labels = "$u(k)$";
end
if nw > 1
    wp_labels = compose("$w_{p,%d}(k)$", 1:nw);
else
    wp_labels = "$w_p(k)$";
end
make_iodplot(Y,Y_m,t,[U Wp],[u_labels wp_labels],{'$y(k), y_m(k)$'},'Time $t$')
filename = "sysid_rodin_stepramp_ioplot";
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

% Plot estimates of selected observer(s)
figure(2); clf
select = find(strcmp("MKF_SP2", obs_labels));
plot_obs_estimates(t, X, Xk_est(:, (select-1)*n+1:select*n), Y, Yk_est(:, select), escape_latex_chars(obs_labels(select)))
filename = "sysid_rodin_stepramp_obs_est_true";
save_fig_to_pdf(fullfile(plot_dir, filename))

% Check results are as expected
assert(all(RMSE_table{:, 'KF3'} ...
    < RMSE_table{:, {'KF1', 'KF2'}}))
assert(all(RMSE_table{:, "MKF_SP1"} ...
    < RMSE_table{:, {'KF1', 'KF2', 'KF3'}}))

% Run experiments to see effect of different disturbance
% model parameters on estimation errors

% Choose simulation case (defined below)
sim_case = 1;

switch sim_case
    case 1
        % Range of paremeter values to test
        epsilon1_values = logspace(-3, -1, 5);
        epsilon2_values = logspace(-3, -1, 5);
        sigma_wp1_values = logspace(-1, 1, 3);
        sigma_wp2_values = logspace(-3, -1, 3);
end

search_space = {
    epsilon1_values, ...
    epsilon2_values, ...
    sigma_wp1_values, ...
    sigma_wp2_values
    };
var_names = {'epsilon1', 'epsilon2', 'sigma_wp_1_2', 'sigma_wp_2_2'};
search_space_idx = cellfun(@(x) 1:length(x), search_space, 'UniformOutput', false);
grid_dims = cellfun(@(x) length(x), search_space, 'UniformOutput', false);

% Search grid
[idx1, idx2, idx3, idx4] = ndgrid(search_space_idx{:});
combs = [reshape(idx1, [], 1) reshape(idx2, [], 1) ...
         reshape(idx3, [], 1) reshape(idx4, [], 1)];
n_combs = size(combs, 1);

% Results filename
filename = sprintf("sysid_rodin_stepramp_rmse_results_SISO_%d.csv", sim_case);

if ~exist(fullfile(results_dir, filename), 'file')

    fprintf("Running %d simulations...\n", n_combs)
    results = array2table(nan(n_combs, 9), ...
        'VariableNames', [{'i1', 'i2', 'i3', 'i4'} var_names {'RMSE'}]);

    % Choose simulations to include in plot
    %selected = true(1, n_combs);
    Xk_est_all = cell(1, n_combs);
    Yk_est_all = cell(1, n_combs);

    for i = 1:n_combs

        eps1_test = search_space{1}(combs(i, 1));
        eps2_test = search_space{2}(combs(i, 2));
        sig12_test = search_space{3}(combs(i, 3));
        sig22_test = search_space{4}(combs(i, 4));

        % Multiple model observer with sequence pruning #1
        label = 'MKF_SP';
        P0 = 1000*eye(n);
        Q0 = diag([q1 0 0]);
        R = sigma_M^2;
        nh = MKF_SP2.nh;  % number of filters
        n_min = MKF_SP2.n_min;  % minimum life of cloned filters
        eps_test = [eps1_test; eps2_test];
        sigma_wp_test = [sigma_wp(1, 1) sig12_test; 
                         sigma_wp(2, 1) sig22_test];
        MKF_SP = MKFObserverSP(model,u_meas,P0,eps_test, ...
            sigma_wp_test,Q0,R,nh,n_min,label);

        % Run simulation
        observers = {MKF_SP};
        f_mkf = 1;

        [Xk_est,Yk_est,DiagPk,MKF_vars] = ...
            run_simulation_obs(Y_m,U,Gamma,[],observers,f_mkf);

        % Add observer estimates to array
        %if selected(i)
        %end
        Xk_est_all{1, i} = Xk_est;
        Yk_est_all{1, i} = Yk_est;

        % Calculate estimation errors
        E_obs = Y - Yk_est;
        rmse = sqrt(sum(E_obs.^2));
        results{i, [{'i1', 'i2', 'i3', 'i4'} var_names {'RMSE'}]} = ...
            [combs(i, :) eps1_test eps2_test sig12_test sig22_test rmse];
        fprintf("Sim %d: (%0.2f, %0.2f, %0.2f, %0.2f) %.3f \n", i, ...
            eps1_test, eps2_test, sig12_test, sig22_test, rmse)

    end

    Xk_est_all = cell2mat(Xk_est_all);
    Yk_est_all = cell2mat(Yk_est_all);

    % Save results
    writetable(results, fullfile(results_dir, filename));
    fprintf("Results saved to file '%s'\n", filename)

    % Make plot of all estimates
    figure(3); clf
    lines = plot(t, Yk_est_all, '-'); hold on
    for i = 1:numel(lines)
        lines(i).Color = [0.1 0.4 1 0.2];
    end
    plot(t, Y, 'k-', 'Linewidth', 2);
    xlabel('Time $t$','Interpreter','latex')
    ylabel('$\hat{y}(k)$','Interpreter','latex')
    grid on

    filename = "sysid_obs_est_grid.png";             
    saveas(gcf, fullfile(plot_dir, filename))

else

    fprintf("Loading simulation results from file '%s'\n", filename)
    results = readtable(fullfile(results_dir, filename));
    assert(size(results, 1) == n_combs);

end

% Find simulation with lowest errors
rmse_min = min(results{:, 'RMSE'});
i_min = find(results{:, 'RMSE'} == rmse_min);

fprintf("Best result: simulation %d of %d\n", i_min, n_combs)
disp(array2table(results{i_min, var_names}, 'VariableNames', var_names))

% Make a grid of RMSE values
rmses = nan(n_epsilon, n_sigma_wp);
for i = 1:n_combs
    comb = results{i, {'x', 'y'}};
    eps_test = results{i, 'epsilon'};
    sig2_test = results{i, 'sigma_wp_2'};
    rmses(comb(1), comb(2)) = results{i, 'RMSE'};
end

% Make contour plot
figure(4); clf
x_pts = log(epsilon_values(combs(:, 1)));
y_pts = log(sigma_wp_values(combs(:, 2)));
plot(x_pts, y_pts, '.'); hold on
set(gca, 'TickLabelInterpreter', 'latex')
X = log(epsilon_values(idx1));
Y = log(sigma_wp_values(idx2));
Z = rmses;
[C,h] = contour(X,Y,Z);
clabel(C,h)
plot(log(epsilon_best), log(sigma_wp_best), 'x', 'LineWidth', 2)
plot(log(epsilon), log(sigma_wp(2)), 'o', 'LineWidth', 2)
text(log(epsilon) + 0.1, log(sigma_wp(2)) - 0.1, sprintf("$%.2f$", rmses(i_min)), 'Interpreter', 'latex')
xlabel("$\log{\epsilon}$", 'Interpreter', 'latex')
ylabel("$\log{\sigma_{w_p,2}}$", 'Interpreter', 'latex')
grid on
legend("simulation results", "contours of RMSE", "lowest RMSE", ...
    "true values", 'Interpreter', 'latex')

filename = sprintf("sysid_rodin_stepramp_contplot_%d", sim_case);
save_fig_to_pdf(fullfile(plot_dir, filename))


% Bayesian optimization
%
% results = bayesopt(fun, vars)
% attempts to find values of vars that minimize fun(vars).