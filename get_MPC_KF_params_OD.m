function [Q, R, Gpred] = get_MPC_KF_params_OD(mpcobj)
% [Q, R, Gpred] = get_MPC_KF_params_OD(mpcobj)
% Determines the parameters of the Kalman filter of the
% MPC object in the case of an output disturbance
% model (which is the default if none specified).
% Also returns the discrete-time prediction model used
% by the estimator.
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

    % Calculate Kalman filter gains
    Gmn = ss(eye(ny),'Ts',Ts);
    B_est = [[Gd.B; zeros(size(God.B,1),size(Gd.B,2))] ... 
             [zeros(size(Gd.B,1),size(God.B,2)); God.B] ...
             [zeros(size(Gd.B,1),size(Gmn.B,2)); 
              zeros(size(God.B,1),size(Gmn,2))]];
    D_est = [Gd.D God.D Gmn.D];
    Q = B_est * B_est';
    R = D_est * D_est';
    N = B_est * D_est';
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