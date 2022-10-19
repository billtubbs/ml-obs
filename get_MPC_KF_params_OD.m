function [Q, R, Gpred] = get_MPC_KF_params_OD(mpcobj)
% [Q, R, Gpred] = get_MPC_KF_params_OD(mpcobj)
% Determines the parameters of the Kalman filter of the
% MPC object in the case of an output disturbance
% model (which is the default if none specified).
% Also returns the discrete-time prediction model used
% by the estimator.
%
% This script is based on this worked-example described 
% in the documentation:
%  - Implement Custom State Estimator Equivalent to Built-In Kalman Filter,
%    https://www.mathworks.com/help/mpc/ug/design-estimator-equivalent-to-mpc-built-in-kf.html
%

    % Get output disturbance model
    God = getoutdist(mpcobj);

    % Make sure there are no input disturbances
    assert(isempty(getindist(mpcobj)))

    % Get plant model
    if isct(mpcobj.Model.Plant)
        % Convert to discrete time model
        Gd = c2d(mpcobj.Model.Plant, mpcobj.Ts);
    else
        Gd = mpcobj.Model.Plant;
    end

    % Determine dimensions of model
    [n, nu, ny, Ts, ~] = check_model(Gd);

    % Augment the plant model with the disturbance model
    A = blkdiag(Gd.A,God.A);
    Bu = [Gd.B; zeros(ny,nu)];
    Cm = [Gd.C God.C];
    D = Gd.D;

    % Prediction model 
    Gpred = ss(A,Bu,Cm,D,Ts);

    % Measurement noise model
    Gmn = ss(eye(ny),'Ts',Ts);

    % Noise input matrix - see equation for w(k)
    B_est = padarray(blkdiag(Gd.B, God.B), [0 size(Gmn.B,2)], 'post');

    % Noise direct transmission matrix - see equation for v(k)
    D_est = [Gd.D God.D Gmn.D];

    % Kalman filter noise covariance parameters
    Q = B_est * B_est';
    R = D_est * D_est';
    N = B_est * D_est';

    % Calculate Kalman filter gains
    G = eye(n+ny);
    H = zeros(ny,n+ny);
    [~, L, ~, M] = kalman(ss(A,[Bu G],Cm,[D H],Ts),Q,R,N);

    % Check gains and model are identical to those of the MPC
    [L1,M1,A1,Cm1,Bu1] = getEstimator(mpcobj);
    assert(max(abs(L - L1), [], [1 2]) < 1e-10)
    assert(max(abs(M - M1), [], [1 2]) < 1e-10)
    assert(max(abs(A - A1), [], [1 2]) < 1e-10)
    assert(max(abs(Cm - Cm1), [], [1 2]) < 1e-10)
    assert(max(abs(Bu1 - Bu), [], [1 2]) < 1e-10)

end