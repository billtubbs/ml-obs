%% Randomly Occurring Deterministic Disturbances
%
% Author: Bill Tubbs, July 2021.
%
% Parameter optimization for MKF filters used in rod_obs_test2.m
%
% Example simulation of various observers on the same system
% as the example in section 5 of MacGregor et al (1984) but with
% the RODD disturbance applied to the input of the process.
%
% Title:  Detection and Estimation of Randomly Occurring 
%         Deterministic Disturbances
% Authors:  Robertson, D. G., Kesavan, P., & Lee, J. H. (1995)
% Publication: Proceedings of 1995 American Control Conference, 
%              ACC'95, 6, 4453?4457.
%
% Dependencies: make sure the following files are in the current path
%   - array2table_with_name.m
%   - array2tablerow.m
%   - combinations.m
%   - combinations_lte.m
%   - get_obs_params.m
%   - kalman_filter.m
%   - kalman_filter_ss.m
%   - kalman_update.m
%   - latex.m
%   - luenberger_filter.m
%   - mkf_filter.m
%   - mkf_filter_RODD.m
%   - n_filters.m
%   - objects2tablerow.m
%   - prob_Gammas.m
%   - prob_gamma.m
%   - rod_obs_rodin_step.m
%   - rod_obs_rodin_step_2x2.m
%   - rod_obs_test2.m
%   - rodd_tf.m
%   - run_simulation_obs.m
%   - sample_random_shocks.m
%   - sys_rodin_step.m
%   - sys_rodin_step_2x2.m
%   - update_KF.m
%   - update_MKF.m
%

clear all; clc

plot_dir = 'plots';
results_dir = 'results';
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end
if ~isfolder(results_dir)
    mkdir(results_dir);
end


%% System model

p_case = 1;
switch p_case
    case 1  % 1st order SISO system with RODD step input disturbance
        p_name = 'rodin_step';

        % Load system and disturbance models from file
        sys_rodin_step

    case 2  % 1st order MIMO system (2x2) with RODD step input disturbances
        p_name = 'rodin_step_2x2sym';
        
        % Load system and disturbance models from file
        sys_rodin_step_2x2sym
        
    otherwise
        error("Value error: p_case")

end

Gd    % Process transfer function Gd(z)
HNd   % Output disturbance ARIMA transfer function HNd(z)
HDd   % Input disturbance RODD transfer function HDd(z)
Gpd   % combined system transfer function (2 inputs, 1 output)
Gpss  % state space model of combined system (2 inputs, 1 output)


% Observers - check observability
Qobs = obsv(A,C);
unobs = length(A) - rank(Qobs);
assert(unobs == 0);


%% Parameter search values

f_values = {3, 5, 10};
m_values = {1, 2, 3};
d_values = {1, 2, 3, 5, 10};

% Random number seed range
seed_values = 1:10;

[f_idx, m_idx, d_idx] = ndgrid(f_values, m_values, d_values);
n_combs = numel(f_idx);

n_filt_max = 1000;
beta_min = 0.8;

for i_comb = 1:n_combs

    %% Define observer

    switch p_case
        case 1  % Multiple model filter

            Q1 = 0.01;
            label = 'MKF';
            P0 = 1000*eye(n);
            Q0 = diag([Q1 1]);
            R = sigma_M^2;
            f = f_idx{i_comb};  % fusion horizon
            m = m_idx{i_comb};  % maximum number of shocks
            d = d_idx{i_comb};  % spacing parameter
            fprintf("Simulation %d : f = %d, m = %d, d = %d\n", i_comb, f, m, d)
            MKF = mkf_filter_RODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
                Q0,R,f,m,d,label);

        case 2  % Multiple model filter 2x2 system

            Q1 = 0.01; Q2 = 0.01;
            label = 'MKF';
            P0 = 1000*eye(n);
            Q0 = diag([Q1 Q2 1 1]);
            R = diag(sigma_M.^2);
            f = f_idx{i_comb};  % fusion horizon
            m = m_idx{i_comb};  % maximum number of shocks
            d = d_idx{i_comb};  % spacing parameter
            fprintf("Simulation %d : f = %d, m = %d, d = %d\n", i_comb, f, m, d)
            MKF = mkf_filter_RODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
                Q0,R,f,m,d,label);

    end

    if MKF.n_filt > n_filt_max
        fprintf("Skipping simulation due to n_filt > %\n", n_filt_max)
        continue
    elseif MKF.beta < beta_min
        fprintf("Skipping simulation due to beta < %.2f\n", beta_min)
        continue
    end


    %% Simulation setup

    % Default simulation settings (may be modified in test cases below)
    t_stop = 1000;
    nT = ceil(t_stop / Ts);
    t = Ts * (0:nT)';


    %% Run repeated simulations of observers and system

    for seed = seed_values

        % Reset random number generator
        rng(seed)

        % before running rod_obs_test_run, you must specify
        % the following variables:
        % observers, test, sim_name
        observers = {MKF};
        sim_name = 'popt';
        test = 1;
        rod_obs_test_run

        % Show selected results
        switch p_case
            case 1
                labels = {'MSE_x1_est_MKF', 'MSE_x2_est_MKF'};
            case 2
                labels = {'MSE_x1_est_MKF', 'MSE_x2_est_MKF', 'MSE_x3_est_MKF', 'MSE_x4_est_MKF'};
        end
        results_table(:, [{'Time', 'seed', 't_stop'} labels])

        % Show plots (turn off for parameter searches)
        % rod_obs_test_plot
    
    end
    
end

disp("Finished")