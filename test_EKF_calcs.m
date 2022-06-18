% Test extended Kalman filter calculations by comparing with
% MATLAB extendedKalmanFilter functions and also with the calculation
% method from GEL-7029 course.
%
% See documentation:
%  - https://www.mathworks.com/help/driving/ug/extended-kalman-filters.html
%

clear all

% Specify path to process model files
addpath('~/process-models/arom3')

% Load parameters from file arom3_params.m
arom3_params

% System dimensions
n = 3;
ny = 2;
nu = 2;
assert(size(x0, 1) == n)
assert(size(p0, 1) == nu)

% Augmented system dimensions
na = n + 2;

% Time step for observer
Ts = 1/60;  % hours

% Observer parameters
P0 = diag([1 25 25 0.5 0.5]);  % see p268
R = diag([1; 100]);  % see p268
Q = diag([0.25; 1; 1; 0.0033; 0.0033]);

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
EKF2 = EKF_observer(na,f,h,params,u_meas,y_meas,dfdx,dhdx,Ts,P0, ...
    Q,R,'EKF2',xa0,y0);
  
% Initialize variables for simulating GEL-7029 EKF
EKF3.Q = Q;
EKF3.R = R;
EKF3.P = P0;
EKF3.xkp1_est = [x0; p0];
EKF3.ykp1_est = arom3_MeasurementFcnRodin2(EKF3.xkp1_est, uk0, Ts, params);

% Measurement data

% Data sample
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

% Load simulation data for testing observers
% filename = 'arom3_sim_benchmark.csv';
% data_dir = 'results';
% sim_data = readtable(fullfile(data_dir,filename), ...
%     'PreserveVariableNames',true);
% Y_m = sim_data{:, {'Y_m_1', 'Y_m_2'}};

nT = size(Y_m, 1) - 1;

% Initialize variables for MATLAB observer calcs
Q = EKF.ProcessNoise;
R = EKF.MeasurementNoise;
Pk = EKF.StateCovariance;
xk_est = EKF.State;

% Do observer calculations
sim_data = nan(nT+1, 2+na*3+3);
k_ind = (0:nT)';
t = Ts*k_ind;
for i = 1:numel(k_ind)
    k = k_ind(i);

    % First measurement y_m(0)
    yk_m = Y_m(i,:)';
    
    % This system has no manipulatable inputs
    uk = [];

    yk_est = arom3_MeasurementFcnRodin2(xk_est,uk,Ts,params);

    % Correct observer states using current measurement
    [CorrectedState, CorrectedStateCovariance] = ...
        correct(EKF, yk_m, uk, Ts, params);

    Hk = arom3_MeasurementJacobianFcnRodin2(xk_est, uk, Ts, params);
    Sk = Hk*Pk*Hk' + R;
    Kk = Pk*Hk'*pinv(Sk);
    xk_est = xk_est + Kk*(yk_m - yk_est);
    Pk = Pk - Kk*Sk*Kk';

    assert(all(abs(xk_est - CorrectedState) < 1e-7, [1 2]))  % TODO: should be < 1e-12
    assert(all(abs(Pk - CorrectedStateCovariance) < 1e-12, [1 2]))

    % Predict state at next sample time
    [PredictedState, PredictedStateCovariance] = predict(EKF, uk, Ts, params);

    Fk = arom3_StateJacobianFcnRodin(xk_est, uk, Ts, params);
    xkp1_est = arom3_StateFcnRodin(xk_est,uk,Ts,params);
    Pkp1 = Fk*Pk*Fk' + Q;

    assert(all(abs(xkp1_est - PredictedState) < 1e-12, [1 2]))  % TODO Something goes wrong here!
    assert(all(abs(Pkp1 - PredictedStateCovariance) < 1e-12, [1 2]))

    % Estimates to be used in next timestep x(k/k-1), P(k/k-1)
    xk_est = xkp1_est;
    Pk = Pkp1;

    % Update my observer
    EKF2 = update_EKF(EKF2, yk_m, uk, Ts);
    assert(all(abs(xkp1_est' - EKF2.xkp1_est') < 1e-12, [1 2]))

    % Update GEL-7029 observer
    EKF3.F = arom3_StateJacobianFcnRodin(EKF3.xkp1_est, uk, Ts, params);
    EKF3.H = arom3_MeasurementJacobianFcnRodin2(EKF3.xkp1_est, uk, Ts, params);
    [EKF3.K, EKF3.P] = ekf_update(EKF3.P,EKF3.F,EKF3.H,EKF3.Q,EKF3.R);
    EKF3.xkp1_est = arom3_StateFcnRodin(EKF3.xkp1_est, uk, Ts, params) + EKF3.K * (yk_m - EKF3.ykp1_est);
    EKF3.ykp1_est = arom3_MeasurementFcnRodin2(EKF3.xkp1_est, uk, Ts, params);

    % Save results
    sim_data(i,:) = [k t(i) xkp1_est' EKF2.xkp1_est' EKF3.xkp1_est' ...
        trace(Pkp1) trace(EKF2.Pkp1) trace(EKF3.P)];

end

obs_labels = {'EKF', 'EKF2', 'EKF3'};
state_labels = {'xkp1_1(k)', 'xkp1_2(k)', 'xkp1_3(k)', 'xkp1_4(k)', 'xkp1_5(k)'};
trP_labels = {'trP'};
var_names = [{'k', 't'} cell(1, na*3) cell(1, 3)];
for j = 1:numel(obs_labels)
    obs_label = obs_labels{j};
    for i = 1:na
        var_names{(j-1)*na+i+2} = sprintf('%s_%s', state_labels{i}, obs_label);
    end
    var_names{na*3+2+j} = sprintf('%s_%s', trP_labels{1}, obs_label);
end
sim_results = array2table(sim_data, 'VariableNames', var_names);

% Check observer estimates match calculated values
assert(all( ...
    abs(sim_results{:,var_names(3:7)} - sim_results{:,var_names(8:12)}) ...
    < 1e-8, [1 2]))
assert(all( ...
    abs(sim_results{:,var_names(3:7)} - sim_results{:,var_names(13:17)}) ...
    < 0.1, [1 2]))
assert(all( ...
    abs(sim_results{:,'trP_EKF'} - sim_results{:,'trP_EKF2'}) ...
    < 1e-8, [1 2]))
assert(all( ...
    abs(sim_results{:,'trP_EKF'} - sim_results{:,'trP_EKF3'}) ...
    < 0.1, [1 2]))

% % Plot state estimates
% plot_obs = {'EKF2', 'EKF3'};
% figure(1); clf
% t_sim = sim_results{:,'t'};
% x_labels = {'$T$','$C_h$','$C_t$', '$k_0$', '$U$'};
% if nT < 100
%     linestyle = 'o-';
% else
%     linestyle = '-';
% end
% for i = 1:na
%     subplot(na,1,i)
%     for j = 1:numel(obs_labels)
%         if ~any(strcmp(obs_labels{j}, plot_obs))
%             continue
%         end
%         label = sprintf('%s_xkp1', obs_labels{j});
%         plot(t_sim, sim_results{:,var_names{(j-1)*na+i+2}}, linestyle); hold on
%     end
%     ylabel(x_labels{i}, 'Interpreter','Latex');
%     grid on
%     title(sprintf('State $x_%d$', i), 'Interpreter','Latex')
%     if i == na
%         xlabel('t (hours)', 'Interpreter','Latex');
%     end
%     legend(plot_obs)
% end
% 
% % Plot trace of covariance matrices
% figure(2); clf
% plot(t_sim, sim_results{:, {'trP_EKF', 'trP_EKF2', 'trP_EKF3'}},'Linewidth',2)
% xlabel('t (hours)', 'Interpreter','Latex');
% ylabel('tr $P(k)$', 'Interpreter','Latex');
% set(gca, 'YScale', 'log')
% grid on
% title('Trace of Covariance Matrix', 'Interpreter','Latex')
% legend(obs_labels)