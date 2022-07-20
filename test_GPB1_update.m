% Test GPB1_update.m

clear all
rng(0)


%% Test calculation - SISO system

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
nj = numel(A);

% Input disturbance variance
%sigma_w = 0.1;
sigma_w = 0;

% Process noise std. dev.
sigma_W = [0; 0];

% Measurement noise std. dev.
sigma_M = 0.1;

% Transition probabilities
epsilon = 0.05;
T = [1-epsilon epsilon; epsilon 1-epsilon];
assert(all(sum(T, 2) == 1))

% Observer parameters (same for all observers)
P0 = 10000;
x0 = 0.5;
Q1 = 0.01;
R1 = 0.1^2;
Q2 = 0.01;
R2 = 0.1^2;
Q = {Q1,Q2};
R = {R1,R2};

% Compare to code from Zhao et al.
Model.A = reshape(cell2mat(A),n,n,nj);
Model.B = reshape(cell2mat(B),n,nu,nj);
Model.C = reshape(cell2mat(C),ny,n,nj);
Model.D = repmat(eye(n),n,n,nj);  % D is used for different purpose
Model.Q = reshape(cell2mat(Q),n,n,nj);
Model.R = reshape(cell2mat(R),ny,ny,nj);
Model.TP = T;
xk_est = x0;
yk = sigma_M .* randn();
u = [0.4; 0.6];
Pk = P0;
% Note: GPB1_estimation model does not include known inputs u(k)
% 'u' here stands for mu which is the posterior likelihood
[x_test,P_test,u_test] = GPB1_estimation(xk_est,Pk,yk,Model,u);

uk = 0;
p_seq_g_Yk = u;

% Prepare cell array of structs to store data
n_filt = nj;
filters = cell(1, n_filt);
for j = 1:n_filt
    filters{j}.xkp1_est = nan(n, 1);
    filters{j}.Pkp1 = nan(n);
    filters{j}.ykp1_est = nan(ny, 1);
end

% Do prediction step first
[Xkp1f_est, Ykp1f_est, Pkp1f] = GPB1_predict(A,B,C,Q,nj, ...
    xk_est,Pk,uk);

% Update step
[xk_est,yk_est,Pk,p_seq_g_Yk] = GPB1_update(A,B,C,Q,R,T,Xkp1f_est, ...
    Ykp1f_est,Pkp1f,yk,p_seq_g_Yk);

% Compare
assert(isequal(x_test, xk_est))
assert(isequal(P_test, Pk))
assert(isequal(u_test, p_seq_g_Yk))


%% Test calculation - 2x2 system

% Sample time
Ts = 1;

% NOTE: this is a previous version of the system with lower
% coupling (-0.2) and epsilon = [0.01; 0.01].

% Discrete time state space model
A = [ 0.8890       0     1 -0.2;
           0  0.8890  -0.2    1;
           0       0     1    0;
           0       0     0    1];
B = [    1 -0.2  0  0;
      -0.2    1  0  0;
         0    0  1  0;
         0    0  0  1];
C = [ 0.1110 0         0  0;
             0  0.1110 0  0];
D = zeros(2, 4);
Gpss = ss(A,B,C,D,Ts);

% Dimensions
n = size(A, 1);
nu = size(B, 2);
ny = size(C, 1);

% Designate measured input and output signals
u_meas = [true; true; false; false];
y_meas = [true; true];

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);
nw = sum(~u_meas);

% Disturbance input
Bw = B(:, ~u_meas);
nw = sum(~u_meas);

% RODD random variable parameters
epsilon = [0.01; 0.01];
sigma_M = [0.1; 0.1];
sigma_wp = [0.01 1; 0.01 1];

% Observer parameters
x0 = [0.5; -0.1; 0.05; -0.05];
Aj = repmat({A}, 1, 3);
Buj = repmat({Bu}, 1, 3);
Cj = repmat({C}, 1, 3);
Duj = repmat({Du}, 1, 3);
nj = numel(Aj);
P0 = 1000*eye(n);
Qj = {diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,1)^2]), ...
      diag([0.01 0.01 sigma_wp(1,2)^2 sigma_wp(2,1)^2]), ...
      diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,2)^2])};
Rj = {diag(sigma_M.^2), diag(sigma_M.^2), diag(sigma_M.^2)};
p_gamma = [1-epsilon epsilon]';
Z = [0 0; 1 0; 0 1];  % combinations
p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
p_gamma = p_gamma ./ sum(p_gamma);  % normalized
T = repmat(p_gamma', 3, 1);

% Compare to code from Zhao et al.
Model.A = reshape(cell2mat(Aj),n,n,nj);
Model.B = reshape(cell2mat(Buj),n,nu,nj);
Model.C = reshape(cell2mat(Cj),ny,n,nj);
Model.D = repmat(eye(n),n,n,nj);  % D is used for different purpose
Model.Q = reshape(cell2mat(Qj),n,n,nj);
Model.R = reshape(cell2mat(Rj),ny,ny,nj);
Model.TP = T;
xk_est = x0;
yk = sigma_M .* randn(ny,1);
u = [0.6; 0.2; 0.2];
Pk = P0;
% Note: GPB1_estimation model does not include known inputs u(k)
% 'u' here stands for mu which is the posterior likelihood
[x_test,P_test,u_test] = GPB1_estimation(xk_est,Pk,yk,Model,u);

uk = zeros(nu,1);
p_seq_g_Yk = u;

% Prepare cell array of structs to store data
n_filt = nj;
filters = cell(1, n_filt);
for j = 1:n_filt
    filters{j}.xkp1_est = nan(n, 1);
    filters{j}.Pkp1 = nan(n);
    filters{j}.ykp1_est = nan(ny, 1);
end

% Do prediction step first
[Xkp1f_est, Ykp1f_est, Pkp1f] = GPB1_predict(Aj,Buj,Cj,Qj,nj, ...
    xk_est,Pk,uk);

% Update step
[xk_est,yk_est,Pk,p_seq_g_Yk] = GPB1_update(Aj,Buj,Cj,Qj,Rj,T,Xkp1f_est, ...
    Ykp1f_est,Pkp1f,yk,p_seq_g_Yk);

% Compare
assert(isequal(x_test, xk_est))
assert(isequal(P_test, Pk))
assert(isequal(u_test, p_seq_g_Yk))