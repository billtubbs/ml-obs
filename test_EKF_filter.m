% Test EKF_filter.m

clear all

addpath('~/process-models/arom3/')

n = 5;
f = @arom3_StateFcnRodin;
h = @arom3_MeasurementFcnRodin2;
u_meas = [false; false];
y_meas = [true; true];
dfdx = @arom3_StateJacobianFcnRodin;
dhdx = @arom3_MeasurementJacobianFcnRodin2;
% see Robertson et al. (1998) for simulation details
Ts = 1;
P0 = [1 25 25 0.5 0.5];
R = [1; 100];
sigma_w = [0.0033; 0.0033];
Q = [0.1; 0.1; 0.1; sigma_w];
label = 'EKF1';
x0 = [742; 463; 537; 5; 6];
obs = EKF_filter(n,f,h,u_meas,y_meas,dfdx,dhdx,Ts,P0,Q,R, ...
    label,x0);

% Check attributes set
assert(obs.n == n)
assert(isequal(obs.f, @arom3_StateFcnRodin))
assert(isequal(obs.h, @arom3_MeasurementFcnRodin2))
assert(isequal(obs.u_meas, u_meas))
assert(isequal(obs.y_meas, y_meas))
assert(isequal(obs.dfdx, @arom3_StateJacobianFcnRodin))
assert(isequal(obs.dhdx, @arom3_MeasurementJacobianFcnRodin2))
assert(obs.Ts == Ts)
assert(isequal(obs.P0, P0))
assert(isequal(obs.Q, Q))
assert(isequal(obs.R, R))
assert(isequal(obs.label, label))
assert(isequal(obs.xkp1_est, x0))

% Test instantiation with unspecified initial condition
obs_0 = EKF_filter(n,f,h,u_meas,y_meas,dfdx,dhdx,Ts,P0,Q,R, ...
    label);
assert(isequal(obs_0.xkp1_est, zeros(5, 1)))



