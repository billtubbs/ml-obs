function [K, P] = ekf_update(P,F,H,Q,R)
% [K, P] = ekf_update(P,F,G,Q,R)
% Extended Kalman filter update
%
    % Compute observer gain (prediction form)
    hphr_inv = pinv(H*P*H' + R);
    K = F*P*H'*hphr_inv;

    % Update covariance matrix
    P = F*(P - P*H'*hphr_inv*H*P)*F' + Q;

end
