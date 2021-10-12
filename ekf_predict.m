function [xkp1, Pkp1] = ekf_predict(xk,Pk,f,Fk,Q,varargin)
% [xkp1, Pkp1] = ekf_predict(xk,Pk,f,Fk,Q,varargin)
% Extended Kalman filter prediction step
%

    % State estimate at next time step
    xkp1 = f(xk, varargin{:});

    % Covariance matrix estimate at next time step
    Pkp1 = Fk*Pk*Fk' + Q;

end