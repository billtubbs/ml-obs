% Test functions luenberger_filter.m and update_KF.m

clear all

% Simulation settings
nT = 200; % number of points
Ts = 3; % sampling period
t = Ts*(0:nT)'; % time vector

% SISO system example from GEL-7029 in file Luenb_no_integr.mlx
A = [0.82 0; 0 0.9];
B = [1; -1];
C = [-0.5 1];
D = 0;
Qp = diag([0.3; 0.2]); % covariance - process noise
Rp = 0.4; % variance - measurement noise

% Dimensions
n = size(A, 1); % number of states
nu = size(B, 2);  % number of inputs
ny = size(C, 1); % number of outputs

% Check if benchmark simulation data file exists
if ~isfile('results/Luenberger_Filter_sim_benchmark.csv')
    error("Run 'Luenb_no_integr_benchmark.mlx' to generate benchmark data.")
end


%% Define and simulate steady-state Kalman filters

% Covariance matrices
Q = diag([0.1; 0.1]);
R = 0.5;
N = zeros(n,ny);

% Define Luenberger observer
poles = [0.9; 0.9];
LB = luenberger_filter(A,B,C,D,Ts,poles,"LB1");

K_test = [0.16; 0];

assert(LB.static_gain == true)
assert(isequal(round(LB.K, 4), K_test))
assert(isequal(LB.xkp1_est, zeros(2, 1)))
assert(LB.ykp1_est == 0)

% seed random number generator
rng(0)

% Measurement noise for the whole simulation
v = sqrt(Rp)*randn(nT,1);

% Process noise for the whole simulation
w = sqrt(Qp)*randn(2,nT);

% Output disturbance
p = zeros(nT,1);
p(t>=300) = 1; % step at time t=300

u0 = 1;  % initial value of u
x0 = inv(eye(length(A)) - A)*B*u0;  % steady-state value of x

% Intialize system (at k = 0)
x = x0;

% Input signal
U = [zeros(10,1); ones(nT+1-10, 1)]; % process input for the whole simulation

% Matrices to collect simulation data
xNprocess = zeros(n, nT+1); % process states
yNprocess = zeros(ny, nT+1); % process outputs
xNobserver = zeros(n, nT+1); % estimated states
yNobserver = zeros(ny, nT+1); % estimated process outputs
xNkalman2 = zeros(n, nT+1); % estimated states
yNkalman2 = zeros(ny, nT+1); % estimated process outputs

for i = 1:nT

    % Process output in current timestep
    y = C*x + v(i) + p(i);
    
    % Record process states and output
    xNprocess(:, i) = x;
    yNprocess(:, i) = y; 

    % Process states in next timestep
    x = A*x + B*U(i) + w(:,i);

    % Lunberger filter update
    LB = update_KF(LB, U(i), y);

    % Record Kalman filter estimates in next timestep
    xNobserver(:, i+1) = LB.xkp1_est;
    yNobserver(:, i+1) = LB.ykp1_est;

end
t = Ts * (0:nT)';


% plot results

% figure(1); clf
% 
% subplot(411);
% plot(t', yNprocess, 'k', t', yNkalman1, 'r', t', yNkalman2, 'g', 'Linewidth', 2)
% legend('Process output', 'KF1 estimates', 'KF2 estimates')
% ylabel('y_1')
% grid on
% 
% subplot(412);
% plot(t', xNprocess(1,:), 'k', t', xNkalman1(1,:), 'r', ...
%     t', xNkalman2(1,:), 'g', 'Linewidth', 2)
% legend('Process state', 'KF1 estimates', 'KF2 estimates')
% ylabel('x_1')
% grid on
% 
% subplot(413);
% plot(t', xNprocess(2,:), 'k', t', xNkalman1(2,:), 'r', ...
%     t', xNkalman2(2,:), 'g', 'Linewidth', 2);
% legend('Process state', 'KF1 estimates', 'KF2 estimates')
% ylabel('x_2')
% grid on
% 
% subplot(414);
% stairs(t', U', 'Linewidth', 2);
% xlabel('Time [s]');
% ylabel('u_1')
% grid on

% Display results
sim_results = [table(t,U) ...
    array2table(xNprocess', 'VariableNames', {'x1', 'x2'}) ...
    array2table(xNobserver', 'VariableNames', {'x1_est', 'x2_est'}) ...
    array2table(yNprocess', 'VariableNames', {'y'}) ...
    array2table(yNobserver', 'VariableNames', {'y_est'}) ...
];

%head(sim_results)

% Verify results by comparing with Luenb_no_integr_benchmark.mlx

filename = 'Luenberger_Filter_sim_benchmark.csv';
bench_sim_results = readtable(fullfile('results', filename));

%head(bench_sim_results)

% Check states
assert(isequal( ...
    round(sim_results{1:200, {'x1', 'x2'}}, 6), ...
    round(bench_sim_results{1:200, {'x1', 'x2'}}, 6) ...
))

% Check state estimates
assert(isequal( ...
    round(sim_results{1:200, {'x1_est', 'x2_est'}}, 6), ...
    round(bench_sim_results{1:200, {'x1_est', 'x2_est'}}, 6) ...
))


% % TODO: Add 2x2 system test
% 
% % Simulate 2x2 system
% 
% % Noise variances
% sigma_p = 0.01;
% sigma_M = 0.1;
% 
% % System model
% A = [0.7 1;
%      0 1];
% B = [1 0;
%      0 1];
% C = [0.3 0];
% D = zeros(1, 2);
% Gpss = ss(A,B,C,D,Ts);
% 
% % Dimensions
% n = size(A, 1);
% nu = size(B, 2);
% ny = size(C, 1);
% 
% %% Define and simulate Kalman filter
% 
% % Discrete time state space model
% Q = sigma_p^2 * diag([1 1]);
% R = sigma_M^2;
% P0 = diag([1e-4 1e-4]);
% KFSS = kalman_filter(A,B,C,D,Ts,P0,Q,R,"KF1");
% 
% % Random inputs
% U = (idinput(size(t)) + 1)/2;
% P = sigma_p*randn(size(t));
% U_sim = [U P];
% 
% % number of points to simulate
% nT = 50;
% t = Ts*(0:nT)';
% 
% % Random inputs
% U = (idinput(size(t)) + 1)/2;
% P = sigma_p*randn(size(t));
% U_sim = [U P];
% 
% [Y, t, X] = lsim(Gpss, U_sim, t);
% 
% 
% X_est = zeros(nT+1, n);
% Y_est = zeros(nT+1, ny);
% E_obs = zeros(nT+1, ny);
% K_obs = zeros(nT+1, n);
% trP_obs = zeros(nT+1, 1);
% 
% xk_est = zeros(n, 1);
% 
% for i=1:nT
%     
%     yk = Y(i,:)';
%     if i > 1
%         uk = U_sim(i-1, 1);
%     else
%         uk = 0;
%     end
%     
%     % Kalman update equations
%     % Update observer gains and covariance matrix
%     KFSS = update_KF(KFSS, [uk; 0], yk);
% 
%     % Record observer estimates
%     X_est(i, :) = KFSS.xkp1_est';
%     Y_est(i, :) = KFSS.ykp1_est';
%     E_obs(i, :) = yk' - KFSS.ykp1_est';
%     K_obs(i, :) = KFSS.K';
%     trP_obs(i, :) = trace(KFSS.P);
% 
% end
