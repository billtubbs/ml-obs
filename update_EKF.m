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
    obs.F = obs.dfdx(obs.xkp1_est, varargin{:});
    obs.H = obs.dhdx(obs.xkp1_est, varargin{:});

    % Previous method according to GEL-7029 prediction form
    [obs.K, obs.P] = ekf_update(obs.P, obs.F, obs.H, obs.Q, obs.R);

    % Update state and output estimates for next timestep
    obs.xkp1_est = obs.f(obs.xkp1_est, varargin{:}) + ...
        obs.K * (yk - obs.ykp1_est);
    obs.ykp1_est = obs.h(obs.xkp1_est, varargin{:});

end