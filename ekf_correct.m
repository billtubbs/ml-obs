function [xk_est, Kk, Pk] = ekf_correct(xk_est,yk_est,yk_m,Pk,Hk,R)
% [xk_est, Kk, Pk] = ekf_correct(xk_est,yk_est,yk_m,Pk,Hk,R)
% Extended Kalman filter correction step
%

    S = Hk*Pk*Hk' + R;

    % Compute observer gain
    Kk = Pk * Hk' * pinv(S);

    % Correct state estimate
    xk_est = xk_est + Kk*(yk_m - yk_est);

    % Update covariance matrix
    % MATLAB docs:
    Pk1 = Pk - Kk*S*Kk';

    % or (see Wikipedia and https://stanford.edu/class/ee363/lectures/ekf.pdf)
    Pk2 = Pk - Kk*Hk*Pk;
    %assert(all(abs(Pk1 - Pk2) < 1e-13, [1 2]))
    %Pk3 = (eye(size(Pk,1)) - Kk*Hk) * Pk;
    %assert(all(abs(Pk2 - Pk3)./(Pk2 + 0.0001*ones(size(Pk2,1))) < 1e-12, [1 2]))
    Pk = Pk2;

end
