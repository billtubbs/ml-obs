function [K, P] = ekf_update(P,F,G,Q,R)
% [K, P] = ekf_update(P,F,G,Q,R)
% Extended Kalman filter update
%
    % Compute observer gain (prediction form)
    hphr_inv = pinv(G*P*G' + R);
    K = F*P*G'*hphr_inv;

    % Update covariance matrix
    P = F*(P - P*G'*hphr_inv*G*P)*F' + Q;

end
