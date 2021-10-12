function [K, P] = ekf_update(P,F,H,Q,R)
% [K, P] = ekf_update(P,F,G,Q,R)
% Extended Kalman filter update - simplified calculation
%
    % Compute observer gain (prediction form)
    S_inv = pinv(H*P*H' + R);
    K = F*P*H'*S_inv;

    % Update covariance matrix
    P = F*(P - P*H'*S_inv*H*P)*F' + Q;

end
