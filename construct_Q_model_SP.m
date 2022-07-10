function [Q, p_gamma] = construct_Q_model_SP(Q0, Bw, alpha, var_wp, nw)
% [Q, p_gamma] = construct_Q_model_SP(Q0, Bw, alpha, var_wp, nw);
% Constructs the parameters needed to model the sub-optimal
% multi-model algorithm using sequence pruning for the
% tracking of infrequently-occurring random disturbances.
% Returns a cell array of the process noise covariance
% matrices Qj for each filter model j = 1:nj for the 
% tracking of infrequently-occurring input disturbances.
%
% Arguments:
%   Q0 : nxn matrix containing variances for measured
%        states on the diagonal (all other elements 
%        are ignored).
%   Bw : system input matrix for the random shock signals
%       (n x nw).
%   alpha : Probability of at least one shock in a 
%       detection interval.
%   var_wp : variances of shock disturbances over 
%       detection interval.
%   nw : number of independent input disturbances.
%

    % TODO: This only works for disturbances with 2
    % states (i.e. shock/no shock). Could be extended
    % to other cases (e.g. no shock, small shock, big
    % shock)
    assert(size(var_wp, 2) == 2)

    % Number of switching models (each with a different Q matrix)
    % This algorithm assumes only one shock per detection period
    % is possible.
    nj = 1 + nw;

    % Array of variances of each shock combination
    var_x = repmat(diag(Q0), 1, nj);

    % Add noise variances for model 1 assuming no shocks
    % occurred
    var_x(:, 1) = var_x(:, 1) + Bw * var_wp(:, 1);

    % Add variances for models 2 to nj to reflect one
    % of each shocks occuring.
    for i = 1:nw
        var_i = var_wp(:, 1);  % no shock
        var_i(i) = var_wp(i, 2);  % shock
        var_x(:, i+1) = var_x(:, i+1) + Bw * var_i;
    end

    Q = cell(1, nj);
    for j = 1:nj
        Q{j} = diag(var_x(:, j));
    end

    % Probabilities of no-shock, shock
    p_gamma = [1-alpha alpha]';

    if nw > 1

        % Possible combinations of each disturbance input:
        % Assume only one may occur in the same sample period
        Z = [zeros(1, nw); eye(nw)];

        % Modified indicator value probabilities
        p_gamma = prod(prob_gamma(Z', p_gamma), 1)';

        % Normalize so that sum(Pr(gamma(k))) = 1
        % TODO: Is this the right thing to do for sub-optimal approach?
        % No I don't think so. If hypotheses don't represent all
        % possible combinations then this is good to know. However,
        % the transition probability matrix must be normalized?
        % p_gamma = p_gamma ./ sum(p_gamma);

    end

end