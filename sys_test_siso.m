% SISO system example from GEL-7029 in file Kalman_Filter.mlx

G = tf(2,conv([10 1],[15 1]));
Ts = 3;
G = ss(c2d(G,Ts,'zoh'));
A = G.a; 
B = G.b;
C = G.c;
D = G.d;

% Dimensions
n = size(A, 1); % number of states
nu = size(B, 2);  % number of inputs
ny = size(C, 1); % number of outputs

% Covariance of process noise
Qp = diag([0.3; 0.2]);

% Variance of measurement noise
Rp = 0.4;

% To check observer convergence use these values
%Qp = [0 0; 0 0]; % process noise
%Rp = 0.0; % measurement noise