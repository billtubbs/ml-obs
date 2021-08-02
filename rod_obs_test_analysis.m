%% Analysis of MKF filter simulations results

% run this file after running rod_obs_test2.m
plot_dir = 'plots';

% Choose MKF observer
f_mkf = 1;
obs = observers{observers_mkf(f_mkf)};

p_seq_g_Yk = sim_out.MKF_p_seq_g_Yk{f_mkf};
% Note: first data points are nans,
% ignore last data point to make plot wider

figure(11); clf
t = Ts*(0:nT)';
ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$Pr(\Gamma(k) \mid Y(k))$'};
filename = sprintf('rod_obs_test2_%d_%s_wfplot.png', test, obs.label);
filepath = fullfile(plot_dir, filename);
show_waterfall_plot(t(2:end-1), p_seq_g_Yk(2:end-1, :), [0 1], ...
    ax_labels, [0 82], filepath);


%% Plot observer shock sequences (Gamma(k))
obs_seq = repmat(cell2mat(obs.seq)', ceil(nT+1 / (obs.d * obs.f)), 1);
obs_seq = obs_seq(1:nT+1, :);

figure(12); clf
p1 = show_waterfall_plot(t(151:251), obs_seq(151:251, :), [0 1], ...
    ax_labels, [0 82]);
hold on
p2 = show_waterfall_plot(t(151:251), p_seq_g_Yk(151:251, :), [0 1], ...
    ax_labels, [0 82], filepath);