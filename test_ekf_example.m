% Test EKF_observer using MATLAB extended Kalman filter examples

clear all

%% Example 1 Van der pol oscilator
%
% From:
%  - https://www.mathworks.com/help/control/ref/extendedkalmanfilter.extendedkalmanfilter.html#bvd_iy8-4

% Initialise MATLAB extendedKalmanFilter object
initialStateGuess = [1; 0];
obj = extendedKalmanFilter(@vdpStateFcn,@vdpMeasurementFcn,initialStateGuess);
obj.StateTransitionJacobianFcn = @vdpStateJacobianFcn;
obj.MeasurementJacobianFcn = @vdpMeasurementJacobianFcn;
obj.ProcessNoise = 0.01;
obj.MeasurementNoise = 0.2;
assert(isequal(obj.StateCovariance, eye(2)))

% Initialize ekf_filter struct
Ts = 1;
P0 = eye(2);
Q = diag([0.01 0.01]);
R = 0.2;
label = 'EKF_vdp';
x0 = [1; 0];
f = @vdpStateFcn;
h = @vdpMeasurementFcn;
u_meas = [];
y_meas = true;
dfdx = @vdpStateJacobianFcn;
dhdx = @vdpMeasurementJacobianFcn;
n = 2; %TODO can get n from P0
obs = EKF_observer(n,f,h,u_meas,y_meas,dfdx,dhdx,Ts,P0,Q,R, ...
    label,x0);

% Check identicial initial parameters
assert(isequal(obs.P, obj.StateCovariance))
assert(isequal(obs.xkp1_est, obj.State))
assert(isequal(obs.Q, obj.ProcessNoise))
assert(isequal(obs.R, obj.MeasurementNoise))
assert(isequal(obj.StateTransitionFcn, obs.f))
assert(isequal(obj.MeasurementFcn, obs.h))
assert(isequal(obj.StateTransitionJacobianFcn, obs.dfdx))
assert(isequal(obj.MeasurementJacobianFcn, obs.dhdx))

% Update extendedKalmanFilter and make prediction
yk_m = 0.9;
[xk2, Pk2_test] = correct(obj, yk_m);
assert(isequal(obj.State, xk2))
assert(isequal(obj.StateCovariance, Pk2_test))
[xkp1_test, Pkp1_test] = predict(obj);
assert(isequal(obj.State, xkp1_test))

% Update ekf_filter and make prediction
obs = update_EKF(obs, yk_m);

% Compare predictions
%assert(isequal(obj.State, obs.xkp1_est))
% 
% ans =
% 
%     0.9167
%    -0.0458
% 
% 
% ans =
% 
%     1.7500
%    -0.0875
%    
% TODO: These are not the same

%assert(isequal(Pkp1_test, obs.P))
% 
% Pkp1_test =
% 
%     0.1792    0.0421
%     0.0421    1.0265
% 
% 
% ans =
% 
%     0.1792    0.0417
%     0.0417    1.0104
%

return




%% Example 2
%
% From:
%  - https://www.mathworks.com/help/control/ref/extendedkalmanfilter.extendedkalmanfilter.html#bvd_iy8-4
% 
% f = @(x,u)(sqrt(x+u));
% h = @(x,v,u)(x+2*u+v^2);