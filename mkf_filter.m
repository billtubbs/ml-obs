function obs = mkf_filter(A,B,C,D,Ts,P0,Q,R,seq,T,label,x0)
% obs = mkf_filter(A,B,C,D,Ts,P0,Q,R,seq,T,label,x0)
%
% Creates a struct for simulating a multi-model Kalman filter
% for state estimation of a Markov jump linear system.
%
% Arguments:
%	A, B, C, D : cell arrays containing discrete-time system
%       matrices for each switching system modelled.
%   Ts : sample period.
%   P0 : cell array of initial covariance matrices of the 
%       state estimates for each filter.
%   Q : cell array of process noise covariance matrices for
%       each switching system.
%   R : cell array of output measurement noise covariance
%       matrices for each switching system.
%   seq : model indicator sequences for each filter (in rows).
%   T : transition probabity matrix of the Markov switching
%       process.
%   label : string name.
%   x0 : intial state estimates (optional).
%

    % Number of switching systems
    nj = numel(A);

    % System dimensions
    [n, nu, ny] = check_dimensions(A{1}, B{1}, C{1}, D{1});
    if nargin == 11
        x0 = zeros(n,1);
    end
    obs.A = A;
    obs.B = B;
    obs.C = C;
    obs.D = D;
    obs.Ts = Ts;
    obs.Q = Q;
    obs.R = R;
    obs.seq = seq;
    obs.T = T;
    obs.label = label;

    % Check transition probability matrix
    assert(all(sum(T, 2) == 1))

    % Initialize covariance matrix of estimation errors
    obs.P = P0;

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
    obs.i = 0;

    % Initialize conditional probabilities: all equally likely
    obs.p_seq_g_Yk = ones(n_filt, 1) ./ n_filt;

    % Empty vectors to store values for filter calculations
    % p(y(k)|Gamma(k),Y(k-1))
    obs.p_yk_g_seq_Ykm1 = zeros(n_filt, 1);
    % Pr(gamma(k)|Y(k-1))
    obs.p_gammak_g_Ykm1 = zeros(n_filt, 1);
    % Pr(gamma(k))
    obs.p_gamma_k = zeros(n_filt, 1);
    % Pr(Gamma(k)|Y(k-1))
    obs.p_seq_g_Ykm1 = zeros(n_filt, 1);

    % Create multi-model observer
    obs.filters = cell(n_filt, 1);
    for i = 1:n_filt
        label_i = sprintf('%s%02d',label,i);
        % Initialize each filter with system #1
        obs.filters{i} = kalman_filter(A{1},B{1},C{1},D{1},Ts,P0{i}, ...
            Q{1},R{1},label_i,x0);
    end

    % Initialize estimates to zero
    obs.xkp1_est = x0;
    obs.ykp1_est = obs.C{1} * obs.xkp1_est;

    % Save useful variables in struct
    obs.nj = nj;
    obs.n = n;
    obs.nu = nu;
    obs.ny = ny;
    obs.n_filt = n_filt;
    obs.nf = nf;

end