% Test calcs made by MATLAB extendedKalmanFilter
%
% See documentation:
%  - https://www.mathworks.com/help/driving/ug/extended-kalman-filters.html
%

clear all

% Specify path to process model files
addpath('~/process-models/arom3')
addpath('~/process-observers')

% Load parameters from file arom3_params.m
arom3_params

% Augmented system dimensions
n = 3;
ny = 2;
nu = 2;
assert(size(x0, 1) == n)
assert(size(p0, 1) == nu)

% Augmented system dimensions
na = n + 2;

% Time step for observer
Ts = 1/60;  % hours

% Augmented model with disturbances as states
P0 = diag([1 25 25 0.5 0.5]);  % see p268
R = diag([1; 100]);  % see p268
cov_w = [0.0033; 0.0033];
% TODO: Doesn't say what process noise on x(1) to x(3) is
Q = diag([0.25; 1; 1; cov_w]);

% EKF observer using MATLAB object
InitialState = [x0; p0];
EKF = extendedKalmanFilter(@arom3_StateFcnRodin, ...
    @arom3_MeasurementFcnRodin2, InitialState);
EKF.MeasurementJacobianFcn = @arom3_MeasurementJacobianFcnRodin2;
EKF.StateTransitionJacobianFcn = @arom3_StateJacobianFcnRodin;
EKF.ProcessNoise = Q;
EKF.MeasurementNoise = R;
EKF.StateCovariance = P0;

% EKF observer using my code
f = @arom3_StateFcnRodin;
h = @arom3_MeasurementFcnRodin2;
u_meas = [false; false];
y_meas = [true; true];
dfdx = @arom3_StateJacobianFcnRodin;
dhdx = @arom3_MeasurementJacobianFcnRodin2;
xa0 = [x0; p0];
uk0 = [];
y0 = arom3_MeasurementFcnRodin2(xa0,uk0,Ts,params);
EKF2 = EKF_observer(na,f,h,u_meas,y_meas,dfdx,dhdx,Ts,P0,Q,R, ...
      'EKF2',xa0,y0);

% This system has no manipulatable inputs
uk = [];

% Measurement data
Y_m = [740.40 521.81
       741.23 536.48
       738.07 539.39
       740.26 541.38
       739.75 524.80
       738.43 537.09
       738.98 533.39
       739.46 541.22
       741.76 538.14
       741.07 542.06];

nT = size(Y_m, 1) - 1;

% Do observer calculations
Q = EKF.ProcessNoise;
R = EKF.MeasurementNoise;
Pk = EKF.StateCovariance;
xk_est = EKF.State;

for k = 0:nT
    i = k + 1;

    % First measurement y_m(0)
    yk_m = Y_m(i,:)';
    yk_est = arom3_MeasurementFcnRodin2(xk_est,uk,Ts,params);

    % Correct observer states using current measurement
    [CorrectedState, CorrectedStateCovariance] = ...
        correct(EKF, yk_m, uk, Ts, params);

    Hk = arom3_MeasurementJacobianFcnRodin2(xk_est, uk, Ts, params);
    Sk = Hk*Pk*Hk' + R;
    Kk = Pk*Hk'*pinv(Sk);
    xk_est = xk_est + Kk*(yk_m - yk_est);
    Pk = Pk - Kk*Sk*Kk';

    assert(all(abs(xk_est - CorrectedState) < 1e-12, [1 2]))
    assert(all(abs(Pk - CorrectedStateCovariance) < 1e-12, [1 2]))

    % Predict state at next sample time
    [PredictedState, PredictedStateCovariance] = predict(EKF, uk, Ts, params);

    Fk = arom3_StateJacobianFcnRodin(xk_est, uk, Ts, params);
    xkp1_est = arom3_StateFcnRodin(xk_est,uk,Ts,params);
    Pkp1 = Fk*Pk*Fk' + Q;

    assert(all(abs(xkp1_est - PredictedState) < 1e-12, [1 2]))
    assert(all(abs(Pkp1 - PredictedStateCovariance) < 1e-12, [1 2]))

    % Estimates to be used in next timestep x(k/k-1), P(k/k-1)
    xk_est = xkp1_est;
    Pk = Pkp1;

    % Update my observer
    EKF2 = update_EKF(EKF2, yk_m, uk, Ts, params);
    assert(all(abs(xkp1_est' - EKF2.xkp1_est') < 1e-12, [1 2]))

end
