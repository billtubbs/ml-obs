function [Kf,xk_est,Pk,yk_est,Sk] = kalman_update_f(C,R,xk_pred,Pk_pred,yk)
% [Kf,xk_est,Pk,yk_est,Sk] = kalman_update_f(C,R,xk_pred,Pk_pred,yk)
% computes the updated correction gain, state estimate, state error 
% covariance, output estimate, and output error covariance in the
% current time instant of a Kalman filter (filtering form) given
% a new measurement:
%
%   S(k) = C(k) * P(k|k-1) * C'(k) + R(k)
%   Kf(k) = P(k)*C(k)'*(C(k)*P(k)*C'(k) + R(k))^-1
%   x_est(k|k) = x_est(k|k-1) + Kf * (y(k) - C(k) * x_est(k|k-1))
%   P(k|k) = P(k) - Kf(k)*S(k)*Kf'(k) + Q(k)
%   y_est(k|k) = C(k) * x_est(k|k)
%
% Arguments:
%   C : (ny, n) matrix
%     Measurement matrix.
%   R : (ny, ny) matrix
%     Measurement error covariance.
%   xk_pred : (n, 1) vector
%     Prior estimate of the states, x_est(k|k-1).
%   Pk_pred : (n, n) matrix
%     Prior estimate of the state error covariance, P(k|k-1).
%   yk : (ny, 1) vector
%     System output measurement at time k.
%
% Returns:
%   Kf : (n, ny) matrix
%     Correction gain matrix.
%   xk_est : (n, 1) vector
%     Posterior estimate of the states, x_est(k|k).
%   Pk : (n, n) matrix
%     State error covariance, P(k|k).
%   yk_est : (ny, 1) vector
%     Posterior estimate of the outputs, y_est(k|k).
%

    % Error covariance of output prediction
    Sk = C * Pk_pred * C' + R;

    % Update correction gain
    Kf = Pk_pred * C' / Sk;

    % Update of prior state estimates using measurements 
    % from current time step to produce 'a posteriori' 
    % state estimates
    xk_est = xk_pred + Kf * (yk - C * xk_pred);

    % Update error covariance of state estimates
    Pk = Pk_pred - Kf * Sk * Kf';

    % Updated output estimate
    yk_est = C * xk_est;

end