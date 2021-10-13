function obs = update_EKF(obs, yk, varargin)
% obs = update_EKF(obs, yk, varargin) updates the gain 
% and covariance matrix of the extended Kalman filter 
% observer and calculates the estimates of the states
% and output at the next sample time.
%
% Arguments:
%   obs : struct containing the Kalman filter variables
%         (see function kalman_filter).
%   yk : vector (ny, 1) of system output measurements at 
%       current sample time k.
%   varargin : provide any additional arguments required
%       by your state transtion function (obs.f).
%
    % Linearize system at current operating point
    % Calculate Jacobian matrices

    % Previous method according to GEL-7029 prediction form
    % which makes both steps
    %[obs.K, obs.P] = ekf_update(obs.P, obs.F, obs.H, obs.Q, obs.R);

    % Add model parameters to arguments cell array
    varargin = [varargin obs.params];

    % Copy estimates calculated in previous timestep
    obs.xk_est = obs.xkp1_est;
    obs.yk_est = obs.ykp1_est;
    obs.P = obs.Pkp1;

    % Correction step - update current estimates based on the new
    % measurement, y(k)
    obs.H = obs.dhdx(obs.xk_est, varargin{:});
    [obs.xk_est, obs.K, obs.P] = ekf_correct(obs.xk_est, obs.yk_est, yk, obs.P, obs.H, obs.R);
    obs.yk_est = obs.h(obs.xk_est, varargin{:});

    % Prediction step - estimate states and covariance matrix in
    % next time step
    obs.F = obs.dfdx(obs.xk_est, varargin{:});
    [obs.xkp1_est, obs.Pkp1] = ekf_predict(obs.xk_est, obs.P, obs.f, obs.F, obs.Q, varargin{:});
    obs.ykp1_est = obs.h(obs.xkp1_est, varargin{:});

end