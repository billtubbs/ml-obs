function MKF = mkf_filter(A,B,C,D,Ts,P0,Q,R,seq,p_gamma,label)
% MKF = mkf_filter(A,B,C,D,Ts,P0,Q,R,seq,p_gamma,label)
%
% Creates a struct for simulating a multi-model Kalman filter
% for state estimation of a Markov jump linear system.
%
% Arguments:
%	A, B, C, D : cell arrays containing discrete time system
%       matrices for each switching system modelled.
%   Ts : sample period.
%   P0 : cell array of initial values of covariance matrices for
%       each filter.
%   Q : cell array of process noise covariance matrices for each
%       switching system.
%   R : cell array of output measurement noise covariance matrices
%       for each switching system.
%   seq : model indicator sequences for each filter (in rows).
%   p_gamma : vector of probabilities of each possible sequence
%       value gamma(k).
%   label : string name.
%

    MKF.A = A;
    MKF.B = B;
    MKF.C = C;
    MKF.D = D;
    MKF.Ts = Ts;
    MKF.Q = Q;
    MKF.R = R;
    MKF.seq = seq;
    MKF.p_gamma = p_gamma;
    MKF.label = label;

    % Initialize covariance matrix of estimation errors
    MKF.P = P0;

    % Number of switching systems
    nj = numel(A);

    % System dimensions
    n = size(A{1}, 1);
    nu = size(B{1}, 2);
    ny = size(C{1}, 1);

    % Check system matrix dimensions. All systems must
    % have same input/output dimensions and number of
    % states.
    for j = 1:nj
        assert(size(A{j}, 2) == n)
        assert(isequal(size(Q{j}), size(A{j})))
        assert(isequal(size(R{j}), [ny ny]))
        assert(size(B{j}, 2) == nu)
        assert(size(C{j}, 1) == ny)
    end

    % Number of filters required
    n_filt = size(seq, 1);

    % Fusion horizon length
    nf = size(cell2mat(seq), 2);

    % Initialize sequence index for online operation
    % Start at 0 because update_MKF increments it before
    % starting the update
    MKF.i = 0;

    % Sequence probabilities Pr(Gamma(k))
    p_seq = prob_Gammas(seq, p_gamma);

    % Tolerance parameter (total probability of defined sequences)
    beta = sum(p_seq);

    % Initialize conditional probabilities Pr(Gamma(k)|Y(k))
    %MKF.p_seq_g_Yk = p_seq ./ beta;
    % Alternative: equal initial probabilities
    MKF.p_seq_g_Yk = ones(n_filt, 1) ./ n_filt;

    % Empty vectors to store values for filter calculations
    % p(y(k)|Gamma(k),Y(k-1))
    MKF.p_yk_g_seq_Ykm1 = zeros(n_filt, 1);

    % Pr(gamma(k)|Y(k-1))
    MKF.p_gammak_g_Ykm1 = zeros(n_filt, 1);

    % Pr(gamma(k))
    obs.p_gamma_k = zeros(n_filt, 1);

    % Pr(Gamma(k)|Y(k-1))
    obs.p_seq_g_Ykm1 = zeros(n_filt, 1);

    % Create multi-model observer
    MKF.filters = cell(n_filt, 1);
    for i = 1:n_filt
        label_i = sprintf('%s%02d',label,i);
        % Initialize each filter with system #1
        MKF.filters{i} = kalman_filter(A{1},B{1},C{1},D{1},Ts,P0{i}, ...
            Q{1},R{1},label_i);
    end

    % Initialize estimates to zero
    MKF.xkp1_est = zeros(n, 1);
    MKF.ykp1_est = zeros(ny, 1);

    % Save useful variables in struct
    MKF.nj = nj;
    MKF.n = n;
    MKF.nu = nu;
    MKF.ny = ny;
    MKF.n_filt = n_filt;
    MKF.nf = nf;
    MKF.beta = beta;
    MKF.p_seq = p_seq;

end