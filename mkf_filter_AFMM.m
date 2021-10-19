function obs = mkf_filter_AFMM(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label,x0)
% obs = mkf_filter_AFMM(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
%     Q0,R,n_filt,f,n_min,label,x0)
%
% Creates a struct for simulating the multi-model Kalman
% filter using the adaptive forgetting through multiple
% models (AFMM) method for state estimation in the presence 
% of infrequently-occurring deterministic disturbances, 
% as described in Eriksson and Isaksson (1996).
%
% Arguments:
%   A, B, C, D : discrete time system matrices.
%   Ts : sample period.
%   u_meas : binary vector indicating measured inputs.
%   P0 : Initial value of covariance matrix of the state
%       estimates.
%   epsilon : probability of a shock disturbance.
%   sigma_wp : standard deviation of shock disturbances.
%   Q0 : initial process noise covariance matrix (n, n) with 
%       variances for each state on the diagonal. The  
%       values Q(i, i) for i representing the unmeasured 
%       input states (where u_meas = false), will
%       be modified for each filter by multiplying by the
%       appropriate variances in sigma_wp.
%   R : output measurement noise covariance matrix (ny, ny).
%   n_filt : number of models (Kalman filters) to utilise.
%   f : length of disturbance sequences to record.
%   n_min : minimum life of cloned filters in number of
%       sample periods.
%   label : string name.
%   x0 : intial state estimates (optional).
%
% Reference:
%  -  Eriksson, P.-G., & Isaksson, A. J. (1996). Classification
%     of Infrequent Disturbances. IFAC Proceedings Volumes, 29(1), 
%     6614–6619. https://doi.org/10.1016/S1474-6670(17)58744-3
%

    % TODO: Could re-introduce spacing parameter d

    % Number of states
    n = check_dimensions(A, B, C, D);

    % Initial state estimates
    if nargin == 15
        x0 = zeros(n,1);
    end

    % Number of input disturbances
    n_dist = sum(~u_meas);
    assert(n_dist > 0);

    % Check there are enough filters in total to accommodate
    % those in the holding group + at least one in main group
    assert(n_min > 0)
    assert(n_filt > 0)
    assert((n_filt - n_dist*n_min) >= n_min, ...
        "ValueError: n_filt is too low.")

    % Process noise covariance matrices for each
    % possible input disturbance
    Q = construct_Q(Q0, B, sigma_wp, u_meas);

    % Number of switching models
    nj = numel(Q);

    % Probabilities of no-shock, shock
    p_gamma = [1-epsilon epsilon]';

    if n_dist > 1

        % Possible combinations of each disturbance input:
        % Assume only one may occur in the same sample period
        Z = [zeros(1, n_dist); eye(n_dist)];
 
        % Modified indicator value probabilities
        p_gamma = prod(prob_gamma(Z', p_gamma), 1)';

        % Normalize so that sum(Pr(gamma(k))) = 1
        % TODO: Is this the right thing to do for sub-optimal approach?
        p_gamma = p_gamma ./ sum(p_gamma);

    end

    % Transition probability matrix
    % Note that for RODD disturbances Pr(gamma(k)) is
    % assumed to be an independent random variable.
    T = repmat(p_gamma', nj, 1);

    % Initialize indicator sequences
    seq = mat2cell(zeros(n_filt, f), ones(1, n_filt), f);

    % Initial index to sequences in main filter group
    % and in holding group
    f_hold = 1:n_dist*n_min;
    f_main = f_hold(end)+1:n_filt;

    % System model doesn't change
    A = repmat({A}, 1, nj);
    B = repmat({B}, 1, nj);
    C = repmat({C}, 1, nj);
    D = repmat({D}, 1, nj);
    R = repmat({R}, 1, nj);

    % Initial covariance matrix is the same for all filters
    P0_init = repmat({P0}, 1, n_filt);

    % Create MKF observer struct
    d = 1;  % TODO: Make this a variable parameter
    obs = mkf_filter(A,B,C,D,Ts,P0_init,Q,R,seq,T,d,label,x0);

    % Add additional variables used by AFMM observer
    obs.f = f;
    obs.n_min = n_min;
    obs.f_main = f_main;
    obs.f_hold = f_hold;
    obs.P0 = P0;
    obs.Q0 = Q0;
    obs.epsilon = epsilon;
    obs.sigma_wp = sigma_wp;
    obs.p_gamma = p_gamma;
    obs.nj = nj;

end


function Q = construct_Q(Q0, B, sigma_wp, u_meas)
% Q = construct_Q(Q0, B, sigma_wp, u_meas) 
% returns a cell array of different process noise 
% covariance matrices Qj for each filter model j = 
% 1:nj for the tracking of infrequently-occurring
% input disturbances.
%
% Arguments:
%   Q0 : nxn matrix containing variances for measured
%        states on the diagonal (all other elements 
%        are ignored).
%   B : system input matrix (n x nu)
%   sigma_wp : standard deviation of shock disturbances.
%   u_meas : binary vector indicating which inputs are
%        measured.
%

    % Number of states
    n = size(B, 1);

    % Check size of initial process covariance matrix
    assert(isequal(size(Q0), [n n]), "ValueError: size(Q0)")

    % Number of inputs
    nu = size(B, 2);
    assert(isequal(size(u_meas), [nu 1]))

    % Number of input disturbances
    n_dist = sum(~u_meas);

    % TODO: This only works for disturbances with 2
    % states (i.e. shock/no shock). Could be extended
    % to other cases (e.g. no shock, small shock, big
    % shock)
    assert(size(sigma_wp, 2) == 2)

    % Number of switching models
    % Assume only one shock at a time possible
    nj = 1 + n_dist;

    % Get noise variances provided for states
    % corresponding to measured inputs.
    var_x = diag(Q0);

    % Set noise variances corresponding to input
    % disturbances to default (no shock) values
    idx = find(~u_meas);
    for i = 1:n_dist
        var_x(idx(i)) = var_x(idx(i)) * ...
            (B(idx(i), idx(i)) * sigma_wp(i,1)).^2;
    end

    % Modify sequences 2:nj for each possible shock
    var_x = repmat(var_x, 1, nj);
    for i = 1:n_dist
        var_x(idx(i), i+1) = (B(idx(i), idx(i)) * sigma_wp(i,2)).^2;
    end

    Q = cell(1, nj);
    for j = 1:nj
        Q{j} = diag(var_x(:,j));
    end

end