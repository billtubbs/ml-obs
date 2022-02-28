% Demonstrate multi-model Kalman filter observer to estimate 
% randomly-occurring deterministic disturbance (RODD).
%

clear all
%show_plots = true;
show_plots = false;

% Generate randomly-occurring shocks
% Reset random number generator
seed = 22;
rng(seed)

% Sequence length
nT = 100;

% RODD random variable parameters
epsilon = 0.01;
sigma_w = [0.01; 1];

[Wp, alpha] = sample_random_shocks(nT+1, epsilon, sigma_w(2), sigma_w(1));

if show_plots
    figure(1)
    subplot(3,1,[1 2])
    stairs(0:nT,Wp); grid on
    ylabel('wp(k)')
    title('Disturbance')
    subplot(3,1,3)
    stairs(0:nT, alpha); grid on
    ylim([-0.1 1.1])
    xlabel('k')
    ylabel('alpha(k)')
end

% Load system model
% First order SISO system with one input disturbance
sys_rodin_step

% The state-space representation of the augmented system from
% above file:
% A, B, C, D, Ts;
% Gpss;

% Save simulation results here
sim_results = struct();

% Simulate system in MATLAB
X0 = zeros(n,1);
t = Ts*(0:nT)';
U = zeros(nT+1,1);
U(t>=5) = 1;
[Y,T,X] = lsim(Gpss,[U Wp],t,X0);
V = sigma_M*randn(nT+1, 1);
Ym = Y + V;  % measurement

if show_plots
    figure(2)
    subplot(2,1,1)
    plot(t,Y,t,Ym); grid on
    ylabel('y(k) and y_m(k)')
    legend('y(k)', 'ym(k)')
    title('Outputs')
    subplot(2,1,2)
    stairs(t, [U Wp]); grid on
    xlabel('k')
    ylabel('u(k) and wp(k)')
    legend('u(k)', 'wp(k)')
    title('Inputs')
end

% Steady-state Kalman filter simulation
% Parameters
Q = diag([0.01^2 0.1^2]);
R = 0.1^2;
Bu = B(:,1);  % observer model without unmeasured inputs
Du = D(:,1);
KFSS = kalman_filter_ss(A,Bu,C,Du,Ts,Q,R,'KFSS');

% Simulate observer
Xk_est = nan(nT+1,n);
Yk_est = nan(nT+1,ny);
obs = KFSS;
for i = 1:nT
    uk = U(i,:)';
    yk = Ym(i,:)';
    obs = update_KF(obs, uk, yk);
    Xk_est(i+1,:) = obs.xkp1_est';
    Yk_est(i+1,:) = obs.ykp1_est';
end

if show_plots
    figure(3)
    subplot(2,1,1)
    plot(t,[Y Yk_est]); grid on
    ylabel('y(k) and y_est(k)')
    legend('y(k)','y_est(k)')
    title('Output')
    subplot(2,1,2)
    plot(t, [X Xk_est]); grid on
    xlabel('k')
    ylabel('xi(k) and xi_est(k)')
    legend('x1(k)','x2(k)','x1_est(k)','x2_est(k)')
    title('States')
end

% Calculate mean-squared error in state estimates
mse_KFSS = mean((X(2:end,:) - Xk_est(2:end,:)).^2, [1 2]);

% Store results
Xk_est_KFSS = Xk_est;
Yk_est_KFSS = Yk_est;

% Kalman filter with time-varying gain
P0 = eye(n);
KF1 = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,'KF1');

% Simulate observer
Xk_est = nan(nT+1,n);
Yk_est = nan(nT+1,ny);
obs = KF1;
for i = 1:nT
    uk = U(i,:)';
    yk = Ym(i,:)';
    obs = update_KF(obs, uk, yk);
    Xk_est(i+1,:) = obs.xkp1_est';
    Yk_est(i+1,:) = obs.ykp1_est';
end

if show_plots
    figure(4)
    subplot(2,1,1)
    plot(t,[Y Yk_est]); grid on
    ylabel('y(k) and y_est(k)')
    legend('y(k)','y_est(k)')
    title('Output')
    subplot(2,1,2)
    plot(t, [X Xk_est]); grid on
    xlabel('k')
    ylabel('xi(k) and xi_est(k)')
    legend('x1(k)','x2(k)','x1_est(k)','x2_est(k)')
    title('States')
end

% Calculate mean-squared error in state estimates
mse_KF1 = mean((X(2:end,:) - Xk_est(2:end,:)).^2, [1 2]);

% Store results
Xk_est_KF1 = Xk_est;
Yk_est_KF1 = Yk_est;

% Sub-optimal multi-model observer 1 simulation
% Sub-optimal multi-model observer as described by Robertson _et al._ (1995).

% Define multi-model filter 1
P0 = eye(n);
Q0 = diag([0.01^2 0]);
R = 0.1^2;
f = 5;  % fusion horizon
m = 1;  % maximum number of shocks
d = 3;  % spacing parameter
MKF1 = mkf_observer_RODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,'MKF1');

% Simulate observer
Xk_est = nan(nT+1,n);
Yk_est = nan(nT+1,ny);
obs = MKF1;
% for debugging
[vec_double, vec_int16] = get_obs_vars_vecs(obs);
dlmwrite(sprintf('test-ml-%s-double.csv', obs.label), vec_double, 'delimiter', ',');
dlmwrite(sprintf('test-ml-%s-int16.csv', obs.label), vec_int16, 'delimiter', ',');    
for i = 1:nT
    uk = U(i,:)';
    yk = Ym(i,:)';

    % For debugging
%     if round(obs.filters{2}.P(2, 2), 4) == 1.3333
%         disp('stop')
%     end
    
    obs = update_MKF(obs, uk, yk);

    % for debugging
    [vec_double, vec_int16] = get_obs_vars_vecs(obs);

    dlmwrite(sprintf('test-ml-%s-double.csv', obs.label), ...
        vec_double, 'delimiter', ',', '-append');
    dlmwrite(sprintf('test-ml-%s-int16.csv', obs.label), ...
        vec_int16, 'delimiter', ',', '-append');

    Xk_est(i+1,:) = obs.xkp1_est';
    Yk_est(i+1,:) = obs.ykp1_est';
end

if show_plots
    figure(5)
    subplot(2,1,1)
    plot(t,[Y Yk_est]); grid on
    ylabel('y(k) and y_est(k)')
    legend('y(k)','y_est(k)')
    title('Output')
    subplot(2,1,2)
    plot(t, [X Xk_est]); grid on
    xlabel('k')
    ylabel('xi(k) and xi_est(k)')
    legend('x1(k)','x2(k)','x1_est(k)',['x2_est(k)'])
    title('States')
end

% Calculate mean-squared error in state estimates
mse_MKF1 = mean((X(2:end,:) - Xk_est(2:end,:)).^2, [1 2]);

% Store results
Xk_est_MKF1 = Xk_est;
Yk_est_MKF1 = Yk_est;


% Sub-optimal multi-model observer 2 simulation
% Sub-optimal multi-model observer as described by Eriksson and Isaksson (1996).

% Define multi-model filter 2
P0 = eye(n);
Q0 = diag([0.01^2 0]);
R = 0.1^2;
f = 10;  % sequence history length
n_filt = 5;  % number of filters
n_min = 2;  % minimum life of cloned filters
MKF2 = mkf_observer_AFMM(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,'MKF2');

% Simulate observer
Xk_est = nan(nT+1,n);
Yk_est = nan(nT+1,ny);
obs = MKF2;
for i = 1:nT
    uk = U(i,:)';
    yk = Ym(i,:)';
    obs = update_AFMM(obs, uk, yk);
    Xk_est(i+1,:) = obs.xkp1_est';
    Yk_est(i+1,:) = obs.ykp1_est';
end

if show_plots
    figure(6)
    subplot(2,1,1)
    plot(t,[Y Yk_est]); grid on
    ylabel('y(k) and y_est(k)')
    legend('y(k)','y_est(k)')
    title('Output')
    subplot(2,1,2)
    plot(t, [X Xk_est]); grid on
    xlabel('k')
    ylabel('xi(k) and xi_est(k)')
    legend('x1(k)','x2(k)','x1_est(k)',['x2_est(k)'])
    title('States')
end

% Calculate mean-squared error in state estimates
mse_MKF2 = mean((X(2:end,:) - Xk_est(2:end,:)).^2, [1 2]);

% Check MSE values
assert(round(mse_KFSS, 4) == 0.0873)
assert(round(mse_KF1, 4) == 0.0917)
assert(round(mse_MKF1, 4) == 0.1029)
assert(round(mse_MKF2, 4) == 0.0614)

% Store results
Xk_est_MKF2 = Xk_est;
Yk_est_MKF2 = Yk_est;

% Combine and save all results to csv file
% These are used by test_MKF_example_sim.m
sim_results = table(t, Xk_est_KFSS, Yk_est_KFSS, Xk_est_KF1, Yk_est_KF1, ...
    Xk_est_MKF1, Yk_est_MKF1, Xk_est_MKF2, Yk_est_MKF2);
filename = 'MKF_example.csv';
results_dir = 'results';
writetable(sim_results, fullfile(results_dir, filename));
