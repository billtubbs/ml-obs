% Test Minimal example used in ReadMe

%TODO: Need to replace KalmanFilter with new KalmanFilterF

clear all

% Known inputs
U = [     0     0     1     1     1     1     1     1     1     1 ...
          1     1     1     1     1     1     1     1     1     1 ...
          1]';

% Output measurements
Ym =  [    0.2688    0.9169   -1.1294    0.7311    0.6694 ...
           0.0032    0.5431    1.0032    2.6715    2.3024 ...
           0.2674    2.4771    1.3345    0.9487    1.3435 ...
           0.8878    0.9311    1.7401    1.7012    1.7063 ...
           1.3341]';

% Sampling period
Ts = 0.5;

% Discrete-time transfer function
Gpd = tf(0.3, [1 -0.7], Ts);

% State-space representation of above process model
A = 0.7;
B = 1;
C = 0.3;
D = 0;

% Kalman filter parameters
P0 = 1;  % estimated variance of the initial state estimate
Q = 0.01;  % estimated process noise variance
R = 0.5^2;  % estimated measurement noise variance
KF1 = KalmanFilter(A,B,C,Ts,P0,Q,R,'KF1');

% Above to be replaced with:

% State-space representation of above process model
model.A = 0.7;
model.B = 1;
model.C = 0.3;
model.Ts = Ts;

% Kalman filter parameters
P0 = 1;  % estimated variance of the initial state estimate
model.Q = 0.01;  % estimated process noise variance
model.R = 0.5^2;  % estimated measurement noise variance
KF2 = KalmanFilterF(model,P0,'KF2');


%% Simulate the observer and record the output estimates:

% Number of sample periods
nT = size(Ym, 1) - 1;
% Array to store observer estimates
Yk_est = nan(nT+1, 1);
Ykp1_est = nan(nT+1, 1);
% Save initial estimate (at t=0)
Ykp1_est(1,:) = KF1.ykp1_est;
for i = 1:nT

    % Update observer with measurements
    KF1.update(Ym(i), U(i));
    KF2.update(Ym(i), U(i));

    % Get estimate of output at next sample time
    Ykp1_est(i+1,:) = KF2.ykp1_est;
    Yk_est(i,:) = KF2.yk_est;

    % Check results identical
    assert(abs(KF1.xkp1_est - KF2.xkp1_est) < 1e-14)
    assert(abs(KF1.Pkp1 - KF2.Pkp1) < 1e-14)

end

% Check results

% Y_est_test = [
%          0    0.1876    0.2997    0.3680    0.5749    0.7048    0.7832 ...
%     0.8460    0.8933    0.9359    0.9625    0.9702    0.9866    0.9924 ...
%     0.9944    0.9978    0.9979    0.9982    1.0024    1.0051    1.0070 ...
% ]';
% assert(isequal(round(Ykp1_est, 4), Y_est_test))

Yk_est_test = [
    0.0712    0.1518    0.0350    0.3370    0.5384    0.6685    0.7658 ...
    0.8374    0.8997    0.9398    0.9529    0.9777    0.9868    0.9905 ...
    0.9958    0.9963    0.9969    1.0030    1.0070    1.0098       NaN
]';
assert(isequaln(round(Yk_est, 4), Yk_est_test))

% Plot observer output estimates to measurement data

% figure(1)
% t = Ts*(0:nT)';
% plot(t,Ym,'o',t,Ykp1_est,'o-',t,Yk_est,'-')
% grid on
% xlabel('Time')
% ylabel('Process output')
% legend('Ym','Ypred','Yest')
% title("Observer estimates compared to process measurements")