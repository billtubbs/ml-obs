%% Estimation of multi-model Kalman filter observer to estimate
% randomly-occurring deterministic disturbance (RODD).
%

clear all


%% Generate randomly-occurring shocks
% Generate a sample sequence of the random variable $w_{p,i} \left(k\right)$ 
% described in Robertson et al. (1995).

% Reset random number generator
seed = 22;
rng(seed)

% Sequence length
nT = 100;

% RODD random variable parameters
epsilon = 0.01;
sigma_w = [0.01; 1];

[Wp, alpha] = sample_random_shocks(nT+1, epsilon, sigma_w(2), sigma_w(1));

figure(1)
subplot(3,1,[1 2])
stairs(0:nT,Wp); grid on
ylabel('wp(k)')
subplot(3,1,3)
stairs(0:nT, alpha); grid on
ylim([-0.1 1.1])
xlabel('k')
ylabel('alpha(k)')

%% Load system model
% 

% First order SISO system with one input disturbance
sys_rodin_step

% The state-space representation of the augmented system from
% above file:
A, B, C, D, Ts
Gpss


%% Simulate system

X0 = zeros(n,1);
t = Ts*(0:nT)';
U = zeros(nT+1,1);
U(t>=5) = 1;
[Y,T,X] = lsim(Gpss,[U Wp],t,X0);
V = sigma_M*randn(nT+1, 1);
Ym = Y + V;  % measurement

figure(2)
subplot(2,1,1)
plot(t,Y,t,Ym); grid on
ylabel('y(k) and y_m(k)')
legend('y(k)', 'ym(k)')
subplot(2,1,2)
stairs(t, [U Wp]); grid on
xlabel('k')
ylabel('u(k) and wp(k)')
legend('u(k)', 'wp(k)')


%% Kalman filter simulation
% The function |kalman_filter| can be used to instantiate a struct variable 
% containing the parameters to simulate a standard Kalman filter:

% Parameters
P0 = [ 0.085704    0.044364
       0.044364    0.036832];  % steady-state values
Q = diag([0.01^2 0.1^2]);
R = 0.1^2;
Bu = B(:,1);  % observer model without unmeasured inputs
Du = D(:,1);
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

figure(3)
subplot(2,1,1)
plot(t,[Y Yk_est]); grid on
ylabel('y(k) and y_est(k)')
legend('y(k)','y_est(k)')
subplot(2,1,2)
plot(t, [X Xk_est]); grid on
xlabel('k')
ylabel('xi(k) and xi_est(k)')
legend('x1(k)','x2(k)','x1_est(k)','x2_est(k)')

% Calculate mean-squared error in state estimates
mse_KF = mean((X(2:end,:) - Xk_est(2:end,:)).^2, [1 2])
assert(round(mse_KF, 4) == 0.0873)


%% Sub-optimal multi-model observer 1 simulation
% Sub-optimal multi-model observer as described by Robertson _et al._ (1995).

% Define multi-model filter
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
for i = 1:nT
    uk = U(i,:)';
    yk = Ym(i,:)';
    obs = update_MKF(obs, uk, yk);
    Xk_est(i+1,:) = obs.xkp1_est';
    Yk_est(i+1,:) = obs.ykp1_est';
end

figure(4)
subplot(2,1,1)
plot(t,[Y Yk_est]); grid on
ylabel('y(k) and y_est(k)')
legend('y(k)','y_est(k)')
subplot(2,1,2)
plot(t, [X Xk_est]); grid on
xlabel('k')
ylabel('xi(k) and xi_est(k)')
legend('x1(k)','x2(k)','x1_est(k)',['x2_est(k)'])

% Calculate mean-squared error in state estimates
mse_MKF1 = mean((X(2:end,:) - Xk_est(2:end,:)).^2, [1 2])
assert(round(mse_MKF1, 4) == 0.1029)


%% Sub-optimal multi-model observer 2 simulation
% Sub-optimal multi-model observer as described by Eriksson and Isaksson (1996).

% Define multi-model filter
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

figure(4)
subplot(2,1,1)
plot(t,[Y Yk_est]); grid on
ylabel('y(k) and y_est(k)')
legend('y(k)','y_est(k)')
subplot(2,1,2)
plot(t, [X Xk_est]); grid on
xlabel('k')
ylabel('xi(k) and xi_est(k)')
legend('x1(k)','x2(k)','x1_est(k)',['x2_est(k)'])

% Calculate mean-squared error in state estimates
mse_MKF2 = mean((X(2:end,:) - Xk_est(2:end,:)).^2, [1 2])
assert(round(mse_MKF2, 4) == 0.0614)


%% Simulate same observers in Simulink

% See this Simulink model file:
model = 'MMKF_example_sim';

% Additional parameters for simulink model
inputs.U = [t U];
inputs.V = [t V];
N = zeros(n,ny);
inputs.Wp = [t Wp];

fprintf("Running Simulink simulation...\n")
sim_out = sim(model, 'ReturnWorkspaceOutputs', 'on');

% Check both Kalman filter state estimates are the same
assert(max(abs(sim_out.X_hat_KF.Data - sim_out.X_hat_KF1.Data), [], [1 2]) < 1e-6)

% Check Kalman filter estimates are close to true system states
assert(mean(abs(sim_out.X_hat_KF.Data - sim_out.X.Data), [1 2]) < 0.5)

disp("Simulations complete")
