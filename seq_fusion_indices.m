function [idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj)
% [idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj)
% Generates three vectors of indices used in the sequence
% fusion sub-optimal multi-model observer algorithm 
% described by Robertson and Lee (1998). These vectors 
% determine the branching, mode transition, and sequence
% merging steps for a given set of sequences defined 
% over the fusion horizon.
%
% For more details, see eqn's 22 and 23 and in Robertson and 
% Lee (1998).
%
% Arguments:
%   seq : (nh, f) integer double
%       Set of nh mode sequences to be modelled over the
%       fustion horizon of length f. At every time
%       instant, past sequences are merged into one of
%       these sequences and others are discarded.
%   nj : integer double
%       number of system modes (2 for a single RODD)
%
% Returns:
%   idx_branch : (nb, 1) integer double
%       This column vector determines how the nh hypotheses
%       will be branched (i.e. split) into nb new hypotheses
%       at the next time instant.
%   idx_modes : (nb, 1) integer double
%       This column vector indicates the mode transitions to
%       extend the nb new hypotheses to the next time instant.
%       Note that modes are zero-based (i.e. 0, 1, 2, ...).
%   idx_merge : (nb, 1) integer double
%       This column vector determines how the nb hypotheses
%       will be merged (i.e. combined) into nh hypotheses
%       at the next time instant.
%
% Example:
% >> nj = 2;  % Number of system modes
% >> f = 3;  % Fusion horizon
% >> m = 1;  % Maximum number of shocks
% >> seq = cell2mat(combinations_lte(f, m))
% 
% seq =
% 
%   4×3 int16 matrix
% 
%    0   0   0
%    1   0   0
%    0   1   0
%    0   0   1
% 
% >> [idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);
% >> [idx_branch idx_modes idx_merge]
% 
% ans =
% 
%      1     0     1
%      1     1     4
%      2     0     1
%      2     1     4
%      3     0     2
%      4     0     3
% 
% References:
%  -  Robertson, D. G., & Lee, J. H. (1998). A method for the
%     estimation of infrequent abrupt changes in nonlinear 
%     systems. Automatica, 34(2), 261–270. 
%     https://doi.org/10.1016/S0005-1098(97)00192-1
%

    % TODO: is nj argument needed?  Why not just take 
    % max value from seq?

    % Number of merged hypotheses modelled
    nh = size(seq, 1);

    % Construct the indices for branching and prediction steps
    % for given sequence model
    idx_branch = reshape(repmat((1:nh), nj, 1), [], 1);
    idx_modes = reshape(repmat((0:nj-1)', 1, nh), [], 1);

    % Hypothesis sequences from time k-f to k:
    seq_kmf_to_k = [seq(idx_branch, :) idx_modes];

    % Drop first (oldest) values in sequences to generate
    % hypothesis sequences from time k-f+1 to k:
    seq_kmfp1_to_k = seq_kmf_to_k(:, 2:end);

    % Drop sequences not defined to be modelled (i.e. not 
    % found in seq)
    seq_to_keep = ismember(seq_kmfp1_to_k, seq, 'rows');
    assert(sum(seq_to_keep) == nh, "ValueError: seq cannot " + ...
        "be used recursively")
    seq_kmfp1_to_k = seq_kmfp1_to_k(seq_to_keep, :);
    idx_branch = idx_branch(seq_to_keep, :);
    idx_modes = idx_modes(seq_to_keep, :);

end