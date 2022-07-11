function [Kf, P] = kalman_update_f(P,A,C,Q,R)
% [Kf, P] = kalman_update_f(P,A,C,Q,R) computes the updated
% gain and covariance matrix of a Kalman filter in
% filtering form:
%
% Kf = P*C'*(C*P*C' + R)^-1;
% P = A*(P - P*C'*(C*P*C' + R)^-1*C*P)*A' + Q;
%

    % Observer gain
    Kf = P * C' / (C * P * C' + R);

    % Covariance matrix
    P = A * (P - Kf * C * P) * A' + Q;

end