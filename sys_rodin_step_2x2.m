%% System model definition
%
% Discrete-time transfer function polynomial models for a
% first-order process with a RODD step disturbance at
% the input to the process.
% 
% Usage: run this script from a main script to define
%  - Gd : Process transfer function Gd(z)
%  - HNd : ARIMA noise disturbance transfer function HNd(z)
%  - HDd : Input disturbance RODD transfer function
%  - Gpd : combined system transfer function (2 inputs, 1 output)
%  - Gpss : state space model of combined system
%
% This script is run by run_obs_sim_spec.m during
% automated simulations.
%

%% Discrete transfer function polynomial models

% Sample time
Ts = 1;

% Process model - continuous time
% This is similar to the reservoir system from GEL-7015, HW #3
s = tf('s');
G11 = -0.7 / (1 + 8.5*s);
G12 = -G11;
G21 = 1.5 / (1 + 16*s);
G22 = G21;
Gc = [G11 G12; G21 G22];
Gd = c2d(Gc,Ts,'zoh');

% ARIMA noise process
thetaN0 = 1;
phiN1 = 0.2;
ThetaN = [0 thetaN0];  % direct transmission
PhiN = [1 -phiN1];
HNd1 = tf(ThetaN, conv(PhiN, [1 -1]), Ts);
HNd2 = tf(ThetaN, conv(PhiN, [1 -1]), Ts);
HNd = [HNd1 0; 0 HNd2];

% RODD step disturbance process
ThetaD = 1;
PhiD = 1;
d = 1;

% RODD ramp disturbance process
% ThetaD = 1;
% PhiD = 1;
% d = 2;

% RODD exponential disturbance process
% ThetaD = 1;
% PhiD = [1 -0.5];
% d = 1;

% Combined disturbance model
HDd1 = rodd_tf(ThetaD, PhiD, d, Ts);
HDd2 = rodd_tf(ThetaD, PhiD, d, Ts);
HDd = [HDd1 0; 0 HDd2];

% Combined system transfer functions
Gpd = [Gd series(Gd, HDd)];

% Combined system - discrete time state space model
% Gpss = minreal(ss(Gpd));
% A = Gpss.A;
% B = Gpss.B;
% C = Gpss.C;
% D = Gpss.D;

% Discrete time state space model
A = [ 0.8890       0  1 -1;
           0  0.9394  1  1;
           0       0  1  0;
           0       0  0  1];
B = [ 1 -1  0  0;
      1  1  0  0;
      0  0  1  0;
      0  0  0  1];
C = [-0.07769  0       0  0;
            0  0.09088 0  0];
D = zeros(2, 4);
Gpss = ss(A,B,C,D,Ts);

% Designate measured input and output signals
u_meas = [true; true; false; false];
y_meas = [true; true];
nu = sum(u_meas);
nw = sum(~u_meas);

% Dimensions
n = size(A, 1);
nu = sum(u_meas);
nw = sum(~u_meas);
ny = size(C, 1);

% Model parameter struct used by observers
model = struct();
model.A = A;
model.B = B;
model.C = C;
model.D = D;
model.Ts = Ts;

% Default initial condition
x0 = zeros(n, 1);

% Parameters for random inputs
% RODD random variable parameters
epsilon = [0.01; 0.01];
sigma_wp = [0.01 1; 0.01 1];

% Process noise standard deviation
sigma_W = zeros(n, 1);

% Measurement noise standard deviation
sigma_M = [0.1; 0.1];

% To test observer with no noise disturbances
% sigma_W = zeros(n, 1);
% sigma_M = zeros(ny, 1);

% Initial state of disturbance process
p0 = zeros(nw, 1);