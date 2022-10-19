function [Q, R, Gpred] = get_MPC_KF_params_ID(mpcobj)
% [Q, R, Gpred] = get_MPC_KF_params_ID(mpcobj)
% Determines the parameters of the Kalman filter of the
% MPC object in the case of an input disturbance
% model. Also returns the discrete-time prediction 
% model used by the estimator.
%

    % Get intput disturbance model
    Gid = getindist(mpcobj);

    % Make sure there are no output disturbances
    assert(isempty(getoutdist(mpcobj)))

    % Get plant model
    if isct(mpcobj.Model.Plant)
        % Convert to discrete time model
        Gd = c2d(mpcobj.Model.Plant, mpcobj.Ts);
    else
        Gd = mpcobj.Model.Plant;
    end

    % Get dimensions of system model
    [n, nu, ny, Ts, ~] = check_model(Gd);
    MVs = Gd.InputGroup.Manipulated;
    nu = numel(MVs);
    UDs = Gd.InputGroup.Unmeasured;
    nw = numel(UDs);

    % Get dimensions of disturbance model
    [n_d, nu_d, ny_d, Ts, ~] = check_model(Gid);

    % Augment the plant model with the disturbance model
    A = [Gd.A Gd.B(:, MVs); zeros(nw,n) Gid.A];
    Bu = [Gd.B(:, MVs) zeros(n,nw); zeros(nw,nu) Gid.B];
    Cm = [Gd.C zeros(ny,nu)];
    D = Gd.D;
    % This should produce the same system as
    % Gpred = series([1 0; 0 Gid], Gd)

    % Prediction model 
    Gpred = ss(A,Bu,Cm,D,Ts);

    % Calculate Kalman filter gains
    Gmn = ss(eye(ny, ny),'Ts',Ts);
    B_est = [[Gd.B;zeros(ny,nu)] [zeros(n,nu);Gid.B] zeros(n+ny,nu)];
    D_est = [Gd.D Gid.D Gmn.D];
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