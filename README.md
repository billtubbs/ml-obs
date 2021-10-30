# process-observers

MATLAB scripts for simulating [process observers](https://en.wikipedia.org/wiki/State_observer) for online state estimation and [sensor fusion](https://en.wikipedia.org/wiki/Sensor_fusion).  

Each observer is produced by a separate MATLAB function but they are implemented as [structure arrays](https://www.mathworks.com/help/matlab/ref/struct.html) ('structs') with many similar attribues that can be passed as arguments to the associated update functions producing similar behaviours.

For example, the following statement creates a struct which contains all the data needed to simulate a discrete-time Kalman filter.

```Matlab
KF1 = kalman_filter(A,B,C,D,Ts,P0,Q,R,'KF1');
```

These structs do not have any sophisticated functionality and cannot be used directly in Simulink.  They are intended to be used for research purposes in hand-coded iterative simulations like this:

```Matlab
% Array to store estimates
y_est = nan(nT+1,1);
% Initial estimate (at t=0)
y_est(1,:) = obs.ykp1_est;
for i = 1:nT

    % Update observer with measurements
    KF1 = update_KF(KF1, u(i), y_m(i));

    % Get estimate of output at next sample time
    y_est(i+1,:) = obs.ykp1_est;

end
```


## Contents

Observers currently included:
- [luenberger_filter.m](luenberger_filter.m) - Luenberger observer (with static correction gain) [[2]](#2).
- [kalman_filter.m](kalman_filter.m) - Kalman filter [[1]](#1).
- [kalman_filter_ss.m](kalman_filter_ss.m) - steady-state Kalman filter (with static correction gain).
- [EKF_observer.m](EKF_observer.m) - extended Kalman filter for non-linear systems.

General-purpose multi-model observers:
- [mkf_observer.m](mkf_observer.m) - multi-model Kalman filter observer.
- [MEKF_observer.m](MEKF_observer.m) - multi-model extended Kalman filter observer.

Specialised multi-model observers:
- [mkf_observer_RODD.m](mkf_observer_RODD.m) and [MEKF_observer_RODD.m](MEKF_observer_RODD.m) - multi-model observers for state estimation in the presence of *randomly-occurring deterministic disturbances* (RODDs) as described in Robertson et al. [[3]](#3).
- [mkf_observer_AFMM.m](mkf_observer_AFMM.m) and [MEKF_observer_AFMM.m](MEKF_observer_AFMM.m) - multi-model observers for state estimation in the presence of *infreuently-occurring disturbances* with the adaptive forgetting through multiple models (AFMM) algorithm as described in Eriksson and Isaksson [[4]](#4).

<p align="center">
    <img src="images/obj_struct_basic.png" alt="Diagram of objects and functions" width="75%">
</p>


## Installation

Clone this repository to your local machine and either add the root to your MATLAB path or work within the main folder.


## Minimal example

Suppose you have some input-output measurement data from a process:
```Matlab
% Measured inputs
u = [     0     0     1     1     1     1     1     1     1     1 ...
          1     1     1     1     1     1     1     1     1     1 ...
          1]';

% Output measurements
y_m = [    0.2688    0.9169   -1.1294    0.7311    0.6694 ...
           0.0032    0.5431    1.0032    2.6715    2.3024 ...
           0.2674    2.4771    1.3345    0.9487    1.3435 ...
           0.8878    0.9311    1.7401    1.7012    1.7063 ...
           1.3341]';

% Sampling period
Ts = 0.5;
```

And, suppose you know the following linear model is a good representation
of the process dynamics:

```Matlab
% Discrete-time transfer function
Gpd = tf(0.3, [1 -0.7], Ts);

% State-space representation of above process model
A = 0.7;
B = 1;
C = 0.3;
D = 0;
```

Define a Kalman filter observer for this process:
```Matlab
% Kalman filter parameters
P0 = 1000;  % estimated variance of the initial state estimate
Q = 0.01;  % estimated process noise variance
R = 0.5^2;  % estimated measurement noise variance
obs = kalman_filter(A,B,C,D,Ts,P0,Q,R,'KF1');
```

Simulate the observer and record the output estimates:
```Matlab
% Number of sample periods
nT = size(y_m,1) - 1;
% Array to store observer estimates
y_est = nan(nT,1);
% Save initial estimate (at t=0)
y_est(1,:) = obs.ykp1_est;
for i = 1:nT

    % update observer
    obs = update_KF(obs, u(i), y_m(i));

    % get estimate of output at next sample time
    y_est(i+1,:) = obs.ykp1_est;

end
```

Compare observer output estimates to measurement data
```Matlab
figure(1)
t = Ts*(0:nT)';
plot(t,y_m,'o',t,y_est,'o-')
grid on
xlabel('Time')
ylabel('Process output')
legend('y_m(k)','y_est(k)')
title("Observer estimates compared to process measurements")
```

<img src='images/siso_kf_example_plot.png' width=600>

## Other examples

See the following LiveScripts for more detailed examples of how to use the functions:

- [kalman_example_SISO.mlx](kalman_example_SISO.mlx) - Kalman filter simulation on a simple single-input, single-output system
- [RODD_code_tutorial.mlx](RODD_code_tutorial.mlx) - Kalman filter and multi-model observer examples on a 2x2 multivariable system

## Help

Most of the functions in this repository are well documented.  Use MATLAB's help command to view the help text:

```Matlab
>> help kalman_filter
```
```text
  obs = kalman_filter(A,B,C,D,Ts,P0,Q,R,label,x0)
  Creates a struct for simulating a time-varying 
  Kalman filter (with time-varying gain).
 
  Arguments:
    A, B, C, D : discrete-time system model matrices.
    Ts : sample period.
    P0 : Initial value of covariance matrix of the state
        estimates.
    Q : Process noise covariance matrix.
    R : Output measurement noise covariance matrix.
    label : string name.
    x0 : intial state estimates (optional).

```

## Testing

A number of unit test scripts are included.  You can run all the tests by running the MATLAB `runtests` command from the root directory.

## References

<a id="1">[1]</a> Kalman, R. E. (1960). A New Approach to Linear Filtering and Prediction Problems. Journal of Basic Engineering. 82: 35–45. https://doi.org/10.1115%2F1.3662552.

<a id="2">[2]</a> Luenberger, D., An Introduction to Observers. IEEE Transactions on Automatic Control 1971, 16 (6), 596–602. https://doi.org/10.1109/TAC.1971.1099826.

<a id="3">[3]</a> Robertson, D. G., Kesavan, P., & Lee, J. H. (1995). Detection and estimation of randomly occurring deterministic disturbances. Proceedings of 1995 American Control Conference - ACC 95, 6, 4453–4457. https://doi.org/10.1109/ACC.1995.532779

<a id="4">[4]</a> Eriksson, P.-G., & Isaksson, A. J. (1996). Classification of Infrequent Disturbances. IFAC Proceedings Volumes, 29(1), 6614–6619. https://doi.org/10.1016/S1474-6670(17)58744-3
