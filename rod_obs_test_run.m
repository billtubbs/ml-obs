%% Run a test simulation and save results to file
% e.g. for parameter optimization

% Before running this script, you must specify the following
% variables: observers, test, sim_name

% For example, this script is used by rod_obs_popt.m

addpath('tblvertcat')

%% Count number of MKF observers

n_obs = numel(observers);
n_obs_mkf = 0;
observers_mkf = double.empty(1, 0);
for i=1:n_obs
    if startsWith(observers{i}.label, "MKF")
        n_obs_mkf = n_obs_mkf + 1;
        observers_mkf(n_obs_mkf) = i;
    end
end


%% Define system input signals

% Array for combined input signals
U = zeros(nT+1, nu);

% Create measured input signals
switch sum(u_meas)
    case 1
        % TODO: Why does adding a step here change the results?
        U(t >= 5, u_meas(1)) = -1;  % add a step
    case 2
        U(t >= 5, u_meas(1)) = -1;  % add a step
        U(t >= 150, u_meas(2)) = 1;  % add a step
end

% Generate randomly-occuring shocks
n_dist = sum(~u_meas);
Wp = zeros(nT+1, n_dist);
alpha = zeros(nT+1, n_dist);
for i = 1:n_dist
    [Wp(:, i), alpha(:, i)] = sample_random_shocks(nT+1, epsilon(i), ...
        sigma_wp(i, 2), sigma_wp(i, 1));
end

U(:, ~u_meas) = Wp;

% Measurement noise signals
% TODO: Does measurement noise start at k=0?
V = zeros(nT+1, ny);
V(t >= 0, y_meas) = sigma_M' .* randn(nT+1, ny);

% Process noise signals
W = zeros(nT+1, n);
W(t >= 0, :) = sigma_W' .* randn(nT+1, n);

% Initial state of process at t=0
x0 = zeros(n, 1);


%% Run simulation

% Compute disturbance signals at input to plant
[Pd, ~] = lsim(HDd, Wp, t);

% Simulate whole system and observers
sim_out = run_simulation_obs(Ts, nT, A, B, C, U, alpha, Pd, V, W, ...
    x0, u_meas, y_meas, observers);

% Updated observers
observers = sim_out.observers;

% Display and save simulation results
sim_out.data
filename = sprintf('rod_obs_%s_%d_%d.csv', sim_name, p_case, test);
writetable(sim_out.data, fullfile(results_dir, filename));

t = sim_out.data{:,'t'};
U = sim_out.data{:,'U'};
alpha = sim_out.data{:,'alpha'};
X = sim_out.data{:,'X'};
X_est = sim_out.data{:,'X_est'};
Y = sim_out.data{:,'Y'};
Y_m = sim_out.data{:,'Y_m'};
Y_est = sim_out.data{:,'Y_est'};
E_obs = sim_out.data{:,'E_obs'};

% Save results from multiple model filters (if used)
sim_out.data_MKF = cell(1, n_obs_mkf);
for f = 1:n_obs_mkf

    label = observers{observers_mkf(f)}.label;
    MKF_sim_results = [sim_out.data(:, {'k', 't'}) ...
        array2table_with_name(sim_out.MKF_i{f}, 'i', '_') ...
        array2table_with_name(sim_out.MKF_p_seq_g_Yk{f}, 'p_seq_g_Yk', '_') ...
        array2table_with_name(sim_out.MKF_X_est{f}, 'X_est', '_') ...
    ];

    filename = sprintf('rod_obs_%s_%d_%d_%s.csv', sim_name, p_case, test, label);
    writetable(MKF_sim_results, fullfile(results_dir, filename));

end


%% Prepare labels for tables and plots

% Make output signal labels
y_labels = cell(1, ny);
for i = 1:ny
    y_labels{i} = sprintf("y_%d(k)", i);
end

% Make input signal labels
u_labels = cell(1, nu);
ui = 1; pi = 1;
for i = 1:nu
    if u_meas(i)
        u_labels{i} = sprintf("u_%d(k)", ui);
        ui = ui + 1;
    else
        u_labels{i} = sprintf("p_%d(k)", pi);
        pi = pi + 1;
    end
end

% Make disturbance signal labels
p_labels = cell(1, ny);
for i = 1:size(Pd, 2)
    p_labels{i} = sprintf("p_%d(k)", i);
end

% Make array of observer labels
obs_labels = cell(1, n_obs);
for i=1:n_obs
    obs_labels{i} = observers{i}.label;
end

% Make array of x_est(k) labels
x_labels = cell(1, n);
x_plot_labels = cell(1, n);
for i = 1:n
    x_labels{i} = sprintf('x%d_est', i);
    x_plot_labels{i} = sprintf('$%s_%d(k)$', '\hat(x)', i);
end

% Make array of labels for each observer 'MSE_x_est(k)'
data_labels = cell(1, n_obs*n);
i = 1;
for j = 1:n
    for f = 1:n_obs
        data_labels{i} = sprintf('MSE_%s_%s', x_labels{j}, obs_labels{f});
        i = i + 1;
    end
end


%% Compute mean-squared errors and save results

% Compute errors in observer state estimates
X_errors = repmat(X, 1, n_obs) - X_est;

% Compute mean-squared errors
% Note: First estimates are at k=1
X_mse = mean(X_errors(2:end, :).^2, 1);

% Display summary table
mse_table = array2table(reshape(X_mse, n, n_obs), ...
    'RowNames', x_labels, ...
    'VariableNames', obs_labels ...
)

% Compute errors in MKF filter estimates (if used)
MKF_X_errors = cell(size(sim_out.MKF_X_est));
MKF_X_mse = cell(1, n_obs_mkf);
for f = 1:n_obs_mkf
    obs = observers{observers_mkf(f)};
    MKF_X_mse{f} = size(sim_out.MKF_X_est{f}, 1);
    % Compute errors in multiple filter state estimates
    % Note: First estimates are at k=1
    MKF_X_errors{f} = repmat(X, 1, obs.n_filt) - sim_out.MKF_X_est{f};
    MKF_X_mse{f} = mean(MKF_X_errors{f}(2:end, :).^2, 1);
end

% Compute errors immediately after shocks
n_shocks = sum(alpha, 1)';
shock_period = 10;
[rows, cols] = ind2sub(size(alpha), find(alpha == 1));
X_mses_after_shocks = cell(1, n_obs);
for f = 1:n_obs
    sum_sq_errors = zeros(shock_period, n);
    counts = zeros(1, n);
    for i = 1:sum(n_shocks)
        errors = X_errors(rows(i):min(end, rows(i)+shock_period-1), (1:n)+(f-1)*n);
        sum_sq_errors(1:size(errors,1),:) = sum_sq_errors(1:size(errors,1),:) + errors.^2;
        counts = counts + size(errors, 1);
    end
    shock_mses = sum_sq_errors ./ counts;
    X_mses_after_shocks{f} = shock_mses;
end


%% Combine all parameters and results and add to summary results file

mse_results = array2table(reshape(mse_table.Variables', 1, n_obs*n), ...
    'VariableNames', data_labels);
sys_model_params = [array2tablerow(A, 'A') array2tablerow(A, 'B') ...
    array2tablerow(A, 'C')];
rv_params = objects2tablerow( ...
    containers.Map({'epsilon', 'sigma_wp', 'sigma_M'}, ...
        {epsilon, sigma_wp, sigma_M}) ...
);
obs_params = cell(1, n_obs);
for f = 1:n_obs
    obs_params{f} = get_obs_params(observers{f});
end
obs_params = horzcat(obs_params{:});
results_table = [ ...
    array2tablerow(datetime(), 'Time') ...
    table(test, seed, t_stop, Ts, nT, nu, ny, n) ...
    sys_model_params ...
    array2tablerow(obs_labels, 'obs') ...
    rv_params ...
    obs_params ...
    mse_results ...
];

% Save to csv file
filename = sprintf('rod_obs_%s_%d_%d_summary.csv', sim_name, p_case, test);
if isfile(fullfile(results_dir, filename))
    % Load existing results and combine
    existing_results = readtable(fullfile(results_dir, filename));
    fprintf("Existing results loaded from file: %s\n", filename)
    results_table = tblvertcat(existing_results, results_table);
end

% Save all results to file
writetable(results_table, fullfile(results_dir, filename));
fprintf("Summary results saved to file: %s\n", filename)