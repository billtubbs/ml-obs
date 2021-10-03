function obs = mkf_filter_RODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label,x0)
% obs = mkf_filter_RODD(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
%   Q0,R,f,m,d,label,x0)
%
% Creates a struct for simulating a multi-model Kalman 
% filter for state estimation in the presence of randomly-
% occurring deterministic disturbances (RODDs) as described
% in Robertson et al. (1995).
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
%        variances for each state on the diagonal. The  
%        values Q(i, i) for i representing the unmeasured 
%        input states (where u_meas = false), will
%        be modified for each filter by multiplying by the
%        appropriate variances in sigma_wp.
%   R : output measurement noise covariance matrix (ny, ny).
%   f : fusion horizon (length of disturbance sequences).
%   m : maximum number of disturbances over fusion horizon.
%   d : spacing parameter in number of sample periods.
%   label : string name.
%   x0 : intial state estimates (optional).
%
% Reference:
%  -  Robertson, D. G., Kesavan, P., & Lee, J. H. (1995). 
%     Detection and estimation of randomly occurring 
%     deterministic disturbances. Proceedings of 1995 American
%     Control Conference - ACC?95, 6, 4453-4457. 
%     https://doi.org/10.1109/ACC.1995.532779
%

    % Number of states
    n = check_dimensions(A, B, C, D);
    
    % Initial state estimates
    if nargin == 15
        x0 = zeros(n,1);
    end

    % Number of input disturbances
    n_dist = sum(~u_meas);

    % Number of filters needed
    n_filt = n_filters(m, f, n_dist);

    % Generate indicator sequences
    S = combinations_lte(f*n_dist, m);

    % Probability of no-shock / shock each period
    % TODO: Is this actually used?
    p = (ones(size(epsilon)) - (ones(size(epsilon)) - epsilon).^d)';
    p = [ones(size(p))-p; p];

    % Probabilities of no-shock, shock
    p_gamma = [1-epsilon epsilon]';

    % Check size of process covariance default matrix
    assert(isequal(size(Q0), [n n]))

    if n_dist == 1

        % Number of Q matrices needed
        nj = 2;

        % Generate required Q matrices
        Q = cell(1, nj);
        for i = 1:nj
            var_x = diag(Q0);
            var_x(~u_meas) = var_x(~u_meas) .* sigma_wp(:, i).^2';
            Q{i} = diag(var_x);
        end

    elseif n_dist > 1
        
        % Note: In the case of more than one input disturbance,
        % there may be multiple combinations of disturbances
        % occuring simultaneously. To simulate these, construct
        % a different Q matrix for each possible combination.
        
        % Reshape sequence data so that each row represents a
        % an input disturbance sequence
        S = cellfun(@(x) reshape(x, n_dist, []), S, ...
            'UniformOutput', false);

        % Find unique combinations of simultaneous shocks
        [Z,~,ic] = unique(cell2mat(S')', 'sorted', 'rows');

        % Number of Q matrices needed
        nj = size(Z, 1);

        % Rearrange as one sequence for each filter and convert
        % back to cell array
        S = reshape((ic - 1)', [], n_filt)'; 
        S = mat2cell(S, ones(n_filt, 1), f);

        % Generate required Q matrices
        Q = cell(1, nj);
        for i = 1:nj
            ind = Z(i, :) + 1;
            var_x = diag(Q0);
            idx = sub2ind(size(sigma_wp), 1:n_dist, ind);
            var_x(~u_meas) = var_x(~u_meas) .* sigma_wp(idx).^2';
            Q{i} = diag(var_x);
        end

        % Modified indicator value probabilities
        p_gamma = prod(prob_gamma(Z', p_gamma), 1)';

        % Normalize so that sum(Pr(gamma(k))) = 1
        % TODO: Is this the right thing to do for sub-optimal approach?
        p_gamma = p_gamma ./ sum(p_gamma);

    else
        
        error("Value error: no unmeasured inputs")

    end

    % Transition probability matrix
    % Note that for RODD disturbances Pr(gamma(k)) is
    % assumed to be an independent random variable.
    T = repmat(p_gamma', nj, 1);

    % Create lengthened indicator sequences by inserting
    % zeros between periods when shocks occur.
    seq = cell(n_filt, 1);
    for i = 1:n_filt
        % Fusion horizon
        f = size(S{i}, 2);
        seq{i} = zeros(size(S{i}, 1), f*d);
        seq{i}(:, (1:f) * d) = S{i};
    end

    % Sequence probabilities Pr(Gamma(k))
    p_seq = prob_Gammas(seq, p_gamma);

    % Tolerance parameter (total probability of defined sequences)
    beta = sum(p_seq);

    % System model does not change
    A = repmat({A}, 1, nj);
    B = repmat({B}, 1, nj);
    C = repmat({C}, 1, nj);
    D = repmat({D}, 1, nj);
    R = repmat({R}, 1, nj);

    % Initial covariance matrix is the same for all filters
    P0_init = repmat({P0}, 1, n_filt);

    % Create MKF observer struct
    obs = mkf_filter(A,B,C,D,Ts,P0_init,Q,R,seq,T,label,x0);

    % Add additional variables used by RODD observer
    obs.S = S;
    obs.P0 = P0;
    obs.Q0 = Q0;
    obs.epsilon = epsilon;
    obs.sigma_wp = sigma_wp;
    obs.f = f;
    obs.m = m;
    obs.d = d;
    obs.p = p;
    obs.beta = beta;
    obs.p_gamma = p_gamma;
    obs.p_seq = p_seq;
    obs.nj = nj;

end